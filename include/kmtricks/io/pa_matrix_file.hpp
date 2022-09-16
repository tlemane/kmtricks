/*****************************************************************************
 *   kmtricks
 *   Authors: T. Lemane
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#pragma once
#include <kmtricks/io/io_common.hpp>
#include <kmtricks/kmer.hpp>
#include <kmtricks/utils.hpp>

namespace km {

class PAMatrixFileHeader : public KmHeader
{
public:
  PAMatrixFileHeader() {};

  void serialize(std::ostream* stream)
  {
    _serialize(stream);
    stream->write(reinterpret_cast<char*>(&matrix_magic), sizeof(matrix_magic));
    stream->write(reinterpret_cast<char*>(&kmer_size), sizeof(kmer_size));
    stream->write(reinterpret_cast<char*>(&kmer_slots), sizeof(kmer_slots));
    stream->write(reinterpret_cast<char*>(&bits), sizeof(bits));
    stream->write(reinterpret_cast<char*>(&bytes), sizeof(bytes));
    stream->write(reinterpret_cast<char*>(&id), sizeof(id));
    stream->write(reinterpret_cast<char*>(&partition), sizeof(partition));
  }

  void deserialize(std::istream* stream)
  {
    _deserialize(stream);
    stream->read(reinterpret_cast<char*>(&matrix_magic), sizeof(matrix_magic));
    stream->read(reinterpret_cast<char*>(&kmer_size), sizeof(kmer_size));
    stream->read(reinterpret_cast<char*>(&kmer_slots), sizeof(kmer_slots));
    stream->read(reinterpret_cast<char*>(&bits), sizeof(bits));
    stream->read(reinterpret_cast<char*>(&bytes), sizeof(bytes));
    stream->read(reinterpret_cast<char*>(&id), sizeof(id));
    stream->read(reinterpret_cast<char*>(&partition), sizeof(partition));
  }

  void sanity_check()
  {
    _sanity_check();
    if (matrix_magic != MAGICS.at(KM_FILE::PAMATRIX))
      throw IOError("Invalid file format.");
  }

public:
  uint64_t matrix_magic {MAGICS.at(KM_FILE::PAMATRIX)};
  uint32_t kmer_size;
  uint32_t kmer_slots;
  uint32_t bits;
  uint32_t bytes;
  uint32_t id;
  uint32_t partition;
};

template<size_t buf_size = 8192>
class PAMatrixWriter : public IFile<PAMatrixFileHeader, std::ostream, buf_size>
{
  using ocstream = lz4_stream::basic_ostream<buf_size>;
public:
  PAMatrixWriter(const std::string& path,
               uint32_t kmer_size,
               uint32_t bits,
               uint32_t id,
               uint32_t partition,
               bool lz4)
    : IFile<PAMatrixFileHeader, std::ostream, buf_size>(path, std::ios::out | std::ios::binary)
  {
    this->m_header.compressed = lz4;
    this->m_header.kmer_size = kmer_size;
    this->m_header.kmer_slots = (kmer_size + 31) / 32;
    this->m_header.bits = bits;
    this->m_header.bytes = NBYTES(bits);
    this->m_header.id = id;
    this->m_header.partition = partition;

    this->m_header.serialize(this->m_first_layer.get());

    this->template set_second_layer<ocstream>(this->m_header.compressed);
  }

  template<size_t MAX_K>
  void write(Kmer<MAX_K>& kmer, std::vector<uint8_t>& vec)
  {
    this->m_second_layer->write(reinterpret_cast<const char*>(kmer.get_data64()),
                                this->m_header.kmer_slots*8);
    this->m_second_layer->write(reinterpret_cast<char*>(vec.data()),
                                vec.size()*sizeof(uint8_t));
  }
};

template<size_t buf_size = 8192>
class PAMatrixReader : public IFile<PAMatrixFileHeader, std::istream, buf_size>
{
  using icstream = lz4_stream::basic_istream<buf_size>;
public:
  PAMatrixReader(const std::string& path)
    : IFile<PAMatrixFileHeader, std::istream, buf_size>(path, std::ios::in | std::ios::binary)
  {
    this->m_header.deserialize(this->m_first_layer.get());
    this->m_header.sanity_check();
    this->template set_second_layer<icstream>(this->m_header.compressed);
  }

  template<size_t MAX_K>
  bool read(Kmer<MAX_K>& kmer, std::vector<uint8_t>& vec)
  {
    this->m_second_layer->read(reinterpret_cast<char*>(kmer.get_data64_unsafe()),
                                this->m_header.kmer_slots*8);
    this->m_second_layer->read(reinterpret_cast<char*>(vec.data()),
                                vec.size()*sizeof(uint8_t));
    if (!this->m_second_layer->gcount())
      return false;
    return true;
  }

  template<size_t MAX_K>
  void write_as_text(std::ostream& stream)
  {
    Kmer<MAX_K> kmer; kmer.set_k(this->m_header.kmer_size);
    std::vector<uint8_t> vec(this->m_header.bytes);
    while (read<MAX_K>(kmer, vec))
    {
      stream << kmer.to_string();
      size_t i = 0;
      for (auto& b : vec)
      {
        for (uint8_t s=0; s<8; s++)
        {
          stream << " " << (((b >> s) & 1) ? '1' : '0');
          if (++i == this->m_header.bits)
            break;
        }
        if (i == this->m_header.bits)
          break;
      }
      stream << "\n";
    }
  }

  template<size_t MAX_K>
  void write_kmers(std::ostream& stream)
  {
    Kmer<MAX_K> kmer; kmer.set_k(this->m_header.kmer_size);
    std::vector<uint8_t> vec(this->m_header.bytes);

    while (read<MAX_K>(kmer, vec))
    {
      stream << kmer.to_string() << '\n';
    }
  }
};

template<size_t buf_size = 8192>
using pmr_t = std::shared_ptr<PAMatrixReader<buf_size>>;

class PAHashMatrixFileHeader : public KmHeader
{
public:
  PAHashMatrixFileHeader() {};

  void serialize(std::ostream* stream)
  {
    _serialize(stream);
    stream->write(reinterpret_cast<char*>(&matrix_magic), sizeof(matrix_magic));
    stream->write(reinterpret_cast<char*>(&bits), sizeof(bits));
    stream->write(reinterpret_cast<char*>(&bytes), sizeof(bytes));
    stream->write(reinterpret_cast<char*>(&id), sizeof(id));
    stream->write(reinterpret_cast<char*>(&partition), sizeof(partition));
  }

  void deserialize(std::istream* stream)
  {
    _deserialize(stream);
    stream->read(reinterpret_cast<char*>(&matrix_magic), sizeof(matrix_magic));
    stream->read(reinterpret_cast<char*>(&bits), sizeof(bits));
    stream->read(reinterpret_cast<char*>(&bytes), sizeof(bytes));
    stream->read(reinterpret_cast<char*>(&id), sizeof(id));
    stream->read(reinterpret_cast<char*>(&partition), sizeof(partition));
  }

  void sanity_check()
  {
    _sanity_check();
    if (matrix_magic != MAGICS.at(KM_FILE::PAMATRIX_HASH))
      throw IOError("Invalid file format.");
  }

public:
  uint64_t matrix_magic {MAGICS.at(KM_FILE::PAMATRIX_HASH)};
  uint32_t bits;
  uint32_t bytes;
  uint32_t id;
  uint32_t partition;
};

template<size_t buf_size = 8192>
class PAHashMatrixWriter : public IFile<PAHashMatrixFileHeader, std::ostream, buf_size>
{
  using ocstream = lz4_stream::basic_ostream<buf_size>;
public:
  PAHashMatrixWriter(const std::string& path,
               uint32_t bits,
               uint32_t id,
               uint32_t partition,
               bool lz4)
    : IFile<PAHashMatrixFileHeader, std::ostream, buf_size>(path, std::ios::out | std::ios::binary)
  {
    this->m_header.compressed = lz4;
    this->m_header.bits = bits;
    this->m_header.bytes = NBYTES(bits);
    this->m_header.id = id;
    this->m_header.partition = partition;

    this->m_header.serialize(this->m_first_layer.get());

    this->template set_second_layer<ocstream>(this->m_header.compressed);
  }

  void write(uint64_t hash, std::vector<uint8_t>& vec)
  {
    this->m_second_layer->write(reinterpret_cast<char*>(&hash), sizeof(hash));
    this->m_second_layer->write(reinterpret_cast<char*>(vec.data()),
                                vec.size()*sizeof(uint8_t));
  }
};

template<size_t buf_size = 8192>
class PAHashMatrixReader : public IFile<PAHashMatrixFileHeader, std::istream, buf_size>
{
  using icstream = lz4_stream::basic_istream<buf_size>;
public:
  PAHashMatrixReader(const std::string& path)
    : IFile<PAHashMatrixFileHeader, std::istream, buf_size>(path, std::ios::in | std::ios::binary)
  {
    this->m_header.deserialize(this->m_first_layer.get());
    this->m_header.sanity_check();
    this->template set_second_layer<icstream>(this->m_header.compressed);
  }

  bool read(uint64_t& hash, std::vector<uint8_t>& vec)
  {
    this->m_second_layer->read(reinterpret_cast<char*>(&hash), sizeof(hash));
    this->m_second_layer->read(reinterpret_cast<char*>(vec.data()),
                                vec.size()*sizeof(uint8_t));
    if (!this->m_second_layer->gcount())
      return false;
    return true;
  }

  void write_as_text(std::ostream& stream)
  {
    uint64_t hash;
    std::vector<uint8_t> vec(this->m_header.bytes);
    while (read(hash, vec))
    {
      stream << std::to_string(hash);
      uint32_t i = 0;
      for (auto& b : vec)
      {
        for (uint8_t s=0; s<8; s++)
        {
          stream << " " << (((b >> s) & 1) ? '1' : '0');
          if (++i == this->m_header.bits)
            break;
        }
        if (i == this->m_header.bits)
          break;
      }
      stream << "\n";
    }
  }
};

template<size_t buf_size = 8192>
using phmr_t = std::shared_ptr<PAHashMatrixReader<buf_size>>;

template<size_t MAX_K>
class PAMatrixFileMerger
{
public:
  struct element
  {
    Kmer<MAX_K> value{0};
    std::vector<uint8_t> count;
    bool is_set {false};
  };

public:
  PAMatrixFileMerger(const std::vector<std::string>& paths, uint32_t kmer_size)
    : m_paths(paths), m_kmer_size(kmer_size)
  {
    init_stream();
    init_state();
  }

  const Kmer<MAX_K>& current() const
  {
    return m_current;
  }

  const std::vector<uint8_t>& bits() const
  {
    return m_counts;
  }

  void init_stream()
  {
    for (auto& path: m_paths)
      m_input_streams.push_back(std::make_shared<PAMatrixReader<8192>>(path));
    m_size = m_paths.size();
  }

  void init_state()
  {
    for (size_t i=0; i<m_size; i++)
    {
      m_elements.push_back(element{});
      m_elements[i].value.set_k(m_kmer_size);
      m_elements[i].count.resize(NBYTES(m_input_streams[i]->infos().bits));

      if (read_next(i))
        m_elements[i].is_set = true;

      if ( (!m_current_set || m_elements[i].value < m_current) && m_elements[i].is_set )
      {
        m_current = m_next = m_elements[i].value;
        m_current_set = true;
      }
    }
    m_counts.resize(m_size, 0);
  }

  bool next()
  {
    m_finish = true;
    m_next_set = false;

    m_current = m_next;
    for (size_t i=0; i<m_size; i++)
    {
      if (m_elements[i].is_set && m_elements[i].value == m_current)
      {
        m_finish = false;
        m_counts = m_elements[i].count;

        if (!read_next(i))
          m_elements[i].is_set = false;
      }

      if (m_elements[i].is_set && (!m_next_set || m_elements[i].value < m_next))
      {
        m_next = m_elements[i].value;
        m_next_set = true;
      }
    }
    return !m_finish;
  }

  void write_as_bin(const std::string& path, bool compressed)
  {
    size_t size = m_input_streams[0]->infos().bits;
    PAMatrixWriter mw(path, m_kmer_size, size, 0, -1, compressed);
    while (next())
    {
      mw.template write<MAX_K>(m_current, m_counts);
    }
  }

  void write_as_text(std::ostream& out)
  {
    int bits = m_input_streams[0]->infos().bits;
    while (next())
    {
      out << m_current.to_string();
      int i = 0;
      for (auto& b : m_counts)
      {
        for (uint8_t s=0; s<8; s++)
        {
          out << " " << ( ((b >> s) & 1) ? '1' : '0');
          if (++i == bits)
            break;
        }
        if (i == bits)
          break;
      }
      out << "\n";
    }
  }

  void write_as_text(const std::string& path)
  {
    std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
    write_as_text(out);
  }

  void write_kmers(std::ostream& out)
  {
    while(next())
    {
      out << m_current.to_string() << "\n";
    }
  }

  void write_kmers(const std::string& path)
  {
    std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
    write_kmers(out);
  }


private:
  bool read_next(size_t i)
  {
    return m_input_streams[i]->template read<MAX_K>(m_elements[i].value,
                                                    m_elements[i].count);
  }

private:
  std::vector<std::string> m_paths;

  std::vector<pmr_t<8192>> m_input_streams;
  std::vector<element> m_elements;

  uint32_t m_size;
  uint32_t m_kmer_size;

  Kmer<MAX_K> m_next;
  Kmer<MAX_K> m_current;
  bool m_next_set {false};
  bool m_current_set {false};
  std::vector<uint8_t> m_counts;

  bool m_finish {false};
};

template<size_t MAX_K>
class PAMatrixFileAggregator
{
public:
  PAMatrixFileAggregator(const std::vector<std::string>& paths, uint32_t kmer_size)
    : m_paths(paths), m_kmer_size(kmer_size)
  {

  }

  void write_as_bin(const std::string& path, bool compressed)
  {
    size_t size = PAMatrixReader<8192>(m_paths[0]).infos().bits;
    PAMatrixWriter<8192> kw(path, m_kmer_size, size, 0, -1, compressed);
    Kmer<MAX_K> k; k.set_k(m_kmer_size);
    std::vector<uint8_t> bits(NBYTES(size));
    for (auto& p : m_paths)
    {
      PAMatrixReader<8192> kr(p);
      while (kr.template read<MAX_K>(k, bits))
        kw.template write<MAX_K>(k, bits);
    }
  }

  void write_as_text(std::ostream& out)
  {
    for (auto& p : m_paths)
    {
      PAMatrixReader<8192> kr(p);
      kr.template write_as_text<MAX_K>(out);
    }
  }

  void write_as_text(const std::string& path)
  {
    std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
    write_as_text(out);
  }

  void write_kmers(std::ostream& out)
  {
    for (auto& p : m_paths)
    {
      PAMatrixReader<8192> kr(p);
      kr.template write_kmers<MAX_K>(out);
    }
  }

  void write_kmers(const std::string& path)
  {
    std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
    write_kmers(out);
  }

private:
  std::vector<std::string> m_paths;
  uint32_t m_kmer_size;
};


class PAHashMatrixFileAggregator
{
public:
  PAHashMatrixFileAggregator(const std::vector<std::string>& paths)
    : m_paths(paths)
  {

  }

  void write_as_bin(const std::string& path, bool compressed)
  {
    size_t size = PAHashMatrixReader<8192>(m_paths[0]).infos().bits;
    PAHashMatrixWriter<8192> kw(
      path, size, 0, -1, compressed);
    uint64_t hash;
    std::vector<uint8_t> bits(NBYTES(size));
    for (auto& p : m_paths)
    {
      PAHashMatrixReader<8192> kr(p);
      while (kr.read(hash, bits))
        kw.write(hash, bits);
    }
  }

  void write_as_text(std::ostream& out)
  {
    for (auto& p : m_paths)
    {
      PAHashMatrixReader<8192> kr(p);
      kr.write_as_text(out);
    }
  }

  void write_as_text(const std::string& path)
  {
    std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
    write_as_text(out);
  }

private:
  std::vector<std::string> m_paths;
  uint32_t m_kmer_size;
};

};
