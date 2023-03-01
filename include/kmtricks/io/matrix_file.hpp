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

class MatrixFileHeader : public KmHeader
{
public:
  MatrixFileHeader() {};

  void serialize(std::ostream* stream)
  {
    _serialize(stream);
    stream->write(reinterpret_cast<char*>(&matrix_magic), sizeof(matrix_magic));
    stream->write(reinterpret_cast<char*>(&kmer_size), sizeof(kmer_size));
    stream->write(reinterpret_cast<char*>(&kmer_slots), sizeof(kmer_slots));
    stream->write(reinterpret_cast<char*>(&count_slots), sizeof(count_slots));
    stream->write(reinterpret_cast<char*>(&nb_counts), sizeof(nb_counts));
    stream->write(reinterpret_cast<char*>(&id), sizeof(id));
    stream->write(reinterpret_cast<char*>(&partition), sizeof(partition));
  }

  void deserialize(std::istream* stream)
  {
    _deserialize(stream);
    stream->read(reinterpret_cast<char*>(&matrix_magic), sizeof(matrix_magic));
    stream->read(reinterpret_cast<char*>(&kmer_size), sizeof(kmer_size));
    stream->read(reinterpret_cast<char*>(&kmer_slots), sizeof(kmer_slots));
    stream->read(reinterpret_cast<char*>(&count_slots), sizeof(count_slots));
    stream->read(reinterpret_cast<char*>(&nb_counts), sizeof(nb_counts));
    stream->read(reinterpret_cast<char*>(&id), sizeof(id));
    stream->read(reinterpret_cast<char*>(&partition), sizeof(partition));
  }

  void deserialize(std::istream* stream, bool kasm)
  {
    if (kasm)
    {
      _deserialize(stream);
      matrix_magic = MAGICS.at(KM_FILE::MATRIX);
      uint64_t dummy;
      stream->read(reinterpret_cast<char*>(&dummy), sizeof(dummy));
      stream->read(reinterpret_cast<char*>(&kmer_size), sizeof(kmer_size));
      stream->read(reinterpret_cast<char*>(&kmer_slots), sizeof(kmer_slots));
      stream->read(reinterpret_cast<char*>(&id), sizeof(id));
      stream->read(reinterpret_cast<char*>(&partition), sizeof(partition));
      stream->read(reinterpret_cast<char*>(&count_slots), sizeof(count_slots));
      nb_counts = 1;
    }
    else
    {
      deserialize(stream);
    }
  }

  void sanity_check()
  {
    _sanity_check();
    if (matrix_magic != MAGICS.at(KM_FILE::MATRIX))
      throw IOError("Invalid file format.");
  }

public:
  uint64_t matrix_magic {MAGICS.at(KM_FILE::MATRIX)};
  uint32_t kmer_size;
  uint32_t kmer_slots;
  uint32_t count_slots;
  uint32_t nb_counts;
  uint32_t id;
  uint32_t partition;
};

template<size_t buf_size = 8192>
class MatrixWriter : public IFile<MatrixFileHeader, std::ostream, buf_size>
{
  using ocstream = lz4_stream::basic_ostream<buf_size>;
public:
  MatrixWriter(const std::string& path,
               uint32_t kmer_size,
               uint32_t count_size,
               uint32_t nb_counts,
               uint32_t id,
               uint32_t partition,
               bool lz4)
    : IFile<MatrixFileHeader, std::ostream, buf_size>(path, std::ios::out | std::ios::binary)
  {
    this->m_header.compressed = lz4;
    this->m_header.kmer_size = kmer_size;
    this->m_header.kmer_slots = (kmer_size + 31) / 32;
    this->m_header.count_slots = count_size;
    this->m_header.nb_counts = nb_counts;
    this->m_header.id = id;
    this->m_header.partition = partition;

    this->m_header.serialize(this->m_first_layer.get());

    this->template set_second_layer<ocstream>(this->m_header.compressed);
  }

  template<size_t MAX_K, size_t MAX_C>
  void write(Kmer<MAX_K>& kmer, std::vector<typename selectC<MAX_C>::type>& counts)
  {
    this->m_second_layer->write(reinterpret_cast<const char*>(kmer.get_data64()),
                                this->m_header.kmer_slots*8);
    this->m_second_layer->write(reinterpret_cast<char*>(counts.data()),
                                counts.size()*(requiredC<MAX_C>::value/8));
  }
};

template<size_t buf_size = 8192>
class MatrixReader : public IFile<MatrixFileHeader, std::istream, buf_size>
{
  using icstream = lz4_stream::basic_istream<buf_size>;
public:
  MatrixReader(const std::string& path, bool kasm = false)
    : IFile<MatrixFileHeader, std::istream, buf_size>(path, std::ios::in | std::ios::binary)
  {
    this->m_header.deserialize(this->m_first_layer.get(), kasm);
    this->m_header.sanity_check();

    this->template set_second_layer<icstream>(this->m_header.compressed);
  }

  template<size_t MAX_K, size_t MAX_C>
  bool read(Kmer<MAX_K>& kmer, std::vector<typename selectC<MAX_C>::type>& counts)
  {
    this->m_second_layer->read(reinterpret_cast<char*>(kmer.get_data64_unsafe()),
                                this->m_header.kmer_slots*8);
    this->m_second_layer->read(reinterpret_cast<char*>(counts.data()),
                                counts.size()*(requiredC<MAX_C>::value/8));
    if (!this->m_second_layer->gcount())
      return false;
    return true;
  }

  template<size_t MAX_K, size_t MAX_C>
  bool read(Kmer<MAX_K>& kmer, std::vector<typename selectC<MAX_C>::type>& counts, std::size_t n)
  {
    this->m_second_layer->read(reinterpret_cast<char*>(kmer.get_data64_unsafe()),
                                this->m_header.kmer_slots*8);
    this->m_second_layer->read(reinterpret_cast<char*>(counts.data()),
                                n*(requiredC<MAX_C>::value/8));
    if (!this->m_second_layer->gcount())
      return false;
    return true;
  }

  template<size_t MAX_K, size_t MAX_C>
  void write_as_text(std::ostream& stream)
  {
    Kmer<MAX_K> kmer; kmer.set_k(this->m_header.kmer_size);
    std::vector<typename selectC<MAX_C>::type> counts(this->m_header.nb_counts);
    while (read<MAX_K, MAX_C>(kmer, counts))
    {
      stream << kmer.to_string();
      for (auto& c : counts)
        stream << " " << std::to_string(c);
      stream << "\n";
    }
  }

  template<size_t MAX_K, size_t MAX_C>
  void write_kmers(std::ostream& stream)
  {
    Kmer<MAX_K> kmer; kmer.set_k(this->m_header.kmer_size);
    std::vector<typename selectC<MAX_C>::type> counts(this->m_header.nb_counts);
    while (read<MAX_K, MAX_C>(kmer, counts))
    {
      stream << kmer.to_string() << '\n';
    }
  }
};

class MatrixHashFileHeader : public KmHeader
{
public:
  MatrixHashFileHeader() {};

  void serialize(std::ostream* stream)
  {
    _serialize(stream);
    stream->write(reinterpret_cast<char*>(&matrix_magic), sizeof(matrix_magic));
    stream->write(reinterpret_cast<char*>(&count_slots), sizeof(count_slots));
    stream->write(reinterpret_cast<char*>(&nb_counts), sizeof(nb_counts));
    stream->write(reinterpret_cast<char*>(&id), sizeof(id));
    stream->write(reinterpret_cast<char*>(&partition), sizeof(partition));
  }

  void deserialize(std::istream* stream)
  {
    _deserialize(stream);
    stream->read(reinterpret_cast<char*>(&matrix_magic), sizeof(matrix_magic));
    stream->read(reinterpret_cast<char*>(&count_slots), sizeof(count_slots));
    stream->read(reinterpret_cast<char*>(&nb_counts), sizeof(nb_counts));
    stream->read(reinterpret_cast<char*>(&id), sizeof(id));
    stream->read(reinterpret_cast<char*>(&partition), sizeof(partition));
  }

  void sanity_check()
  {
    _sanity_check();
    if (matrix_magic != MAGICS.at(KM_FILE::MATRIX_HASH))
      throw IOError("Invalid file format.");
  }

public:
  uint64_t matrix_magic {MAGICS.at(KM_FILE::MATRIX_HASH)};
  uint32_t count_slots;
  uint32_t nb_counts;
  uint32_t id;
  uint32_t partition;
};

template<size_t buf_size = 8192>
class MatrixHashWriter : public IFile<MatrixHashFileHeader, std::ostream, buf_size>
{
  using ocstream = lz4_stream::basic_ostream<buf_size>;
public:
  MatrixHashWriter(const std::string& path,
                   uint32_t count_size,
                   uint32_t nb_counts,
                   uint32_t id,
                   uint32_t partition,
                   bool lz4)
    : IFile<MatrixHashFileHeader, std::ostream, buf_size>(path, std::ios::out | std::ios::binary)
  {
    this->m_header.compressed = lz4;
    this->m_header.count_slots = count_size;
    this->m_header.nb_counts = nb_counts;
    this->m_header.id = id;
    this->m_header.partition = partition;

    this->m_header.serialize(this->m_first_layer.get());

    this->template set_second_layer<ocstream>(this->m_header.compressed);
  }

  template<size_t MAX_C>
  void write(uint64_t hash, std::vector<typename selectC<MAX_C>::type>& counts)
  {
    this->m_second_layer->write(reinterpret_cast<char*>(&hash), sizeof(hash));
    this->m_second_layer->write(reinterpret_cast<char*>(counts.data()),
                                counts.size()*(requiredC<MAX_C>::value/8));
  }
};

template<size_t buf_size = 8192>
class MatrixHashReader : public IFile<MatrixHashFileHeader, std::istream, buf_size>
{
  using icstream = lz4_stream::basic_istream<buf_size>;
public:
  MatrixHashReader(const std::string& path)
    : IFile<MatrixHashFileHeader, std::istream, buf_size>(path, std::ios::in | std::ios::binary)
  {
    this->m_header.deserialize(this->m_first_layer.get());
    this->m_header.sanity_check();

    this->template set_second_layer<icstream>(this->m_header.compressed);
  }

  template<size_t MAX_C>
  bool read(uint64_t& hash, std::vector<typename selectC<MAX_C>::type>& counts)
  {
    this->m_second_layer->read(reinterpret_cast<char*>(&hash), sizeof(hash));
    this->m_second_layer->read(reinterpret_cast<char*>(counts.data()),
                                counts.size()*(requiredC<MAX_C>::value/8));
    if (!this->m_second_layer->gcount())
      return false;
    return true;
  }

  template<size_t MAX_C>
  void write_as_text(std::ostream& stream)
  {
    uint64_t hash;
    std::vector<typename selectC<MAX_C>::type> counts(this->m_header.nb_counts);
    while (read<MAX_C>(hash, counts))
    {
      stream << std::to_string(hash);
      for (auto& c : counts)
        stream << " " << std::to_string(c);
      stream << "\n";
    }
  }
};

template<size_t buf_size = 8192>
using mr_t = std::shared_ptr<MatrixReader<buf_size>>;
template<size_t buf_size = 8192>
using mhr_t = std::shared_ptr<MatrixHashReader<buf_size>>;

template<size_t MAX_K, size_t MAX_C>
class MatrixFileMerger
{
  using count_type = typename selectC<MAX_C>::type;
public:
  struct element
  {
    Kmer<MAX_K> value{0};
    std::vector<count_type> count;
    bool is_set {false};
  };

public:
  MatrixFileMerger(const std::vector<std::string>& paths, uint32_t kmer_size)
    : m_paths(paths), m_kmer_size(kmer_size)
  {
    init_stream();
    init_state();
  }

  const Kmer<MAX_K>& current() const
  {
    return m_current;
  }

  const std::vector<count_type>& counts() const
  {
    return m_counts;
  }

  void init_stream()
  {
    for (auto& path: m_paths)
      m_input_streams.push_back(std::make_shared<MatrixReader<8192>>(path));
    m_size = m_paths.size();
  }

  void init_state()
  {
    for (size_t i=0; i<m_size; i++)
    {
      m_elements.push_back(element{});
      m_elements[i].value.set_k(m_kmer_size);
      m_elements[i].count.resize(m_input_streams[i]->infos().nb_counts);

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
    size_t size = m_input_streams[0]->infos().nb_counts;
    MatrixWriter mw(path, m_kmer_size, requiredC<MAX_C>::value/8, size, 0, -1, compressed);
    while (next())
    {
      mw.template write<MAX_K, MAX_C>(m_current, m_counts);
    }
  }

  void write_as_text(std::ostream& out)
  {
    while (next())
    {
      out << m_current.to_string();
      for (auto& c: m_counts)
        out << " " << std::to_string(c);
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
    while (next())
    {
      out << m_current.to_string() << '\n';
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
    return m_input_streams[i]->template read<MAX_K, MAX_C>(m_elements[i].value,
                                                           m_elements[i].count);
  }

private:
  std::vector<std::string> m_paths;

  std::vector<mr_t<8192>> m_input_streams;
  std::vector<element> m_elements;

  uint32_t m_size;
  uint32_t m_kmer_size;

  Kmer<MAX_K> m_next;
  Kmer<MAX_K> m_current;
  bool m_next_set {false};
  bool m_current_set {false};
  std::vector<count_type> m_counts;

  bool m_finish {false};
};

template<size_t MAX_K, size_t MAX_C>
class MatrixFileAggregator
{
public:
  MatrixFileAggregator(const std::vector<std::string>& paths, uint32_t kmer_size)
    : m_paths(paths), m_kmer_size(kmer_size)
  {

  }

  void write_as_bin(const std::string& path, bool compressed)
  {
    size_t size = MatrixReader<8192>(m_paths[0]).infos().nb_counts;
    MatrixWriter<8192> kw(path, m_kmer_size, requiredC<MAX_C>::value/8, size, 0, -1, compressed);
    Kmer<MAX_K> k; k.set_k(m_kmer_size);
    std::vector<typename selectC<MAX_C>::type> counts(size);
    for (auto& p : m_paths)
    {
      MatrixReader<8192> kr(p);
      while (kr.template read<MAX_K, MAX_C>(k, counts))
        kw.template write<MAX_K, MAX_C>(k, counts);
    }
  }

  void write_as_text(std::ostream& out)
  {
    for (auto& p : m_paths)
    {
      MatrixReader<8192> kr(p);
      kr.template write_as_text<MAX_K, MAX_C>(out);
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
      MatrixReader<8192> kr(p);
      kr.template write_kmers<MAX_K, MAX_C>(out);
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


template<size_t MAX_C>
class MatrixHashFileAggregator
{
public:
  MatrixHashFileAggregator(const std::vector<std::string>& paths)
    : m_paths(paths)
  {

  }

  void write_as_bin(const std::string& path, bool compressed)
  {
    size_t size = MatrixHashReader<8192>(m_paths[0]).infos().nb_counts;
    MatrixHashWriter<8192> kw(path, requiredC<MAX_C>::value/8, size, 0, -1, compressed);
    uint64_t hash;
    std::vector<typename selectC<MAX_C>::type> counts(size);
    for (auto& p : m_paths)
    {
      MatrixHashReader<8192> kr(p);
      while (kr.template read<MAX_C>(hash, counts))
        kw.template write<MAX_C>(hash, counts);
    }
  }

  void write_as_text(std::ostream& out)
  {
    for (auto& p : m_paths)
    {
      MatrixHashReader<8192> kr(p);
      kr.template write_as_text<MAX_C>(out);
    }
  }

  void write_as_text(const std::string& path)
  {
    std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
    write_as_text(out);
  }
private:
  std::vector<std::string> m_paths;
};

};

