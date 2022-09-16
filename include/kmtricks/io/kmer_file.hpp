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

class KmerFileHeader : public KmHeader
{
public:
  KmerFileHeader() {};

  void serialize(std::ostream* stream) override
  {
    _serialize(stream);
    stream->write(reinterpret_cast<char*>(&kmer_magic), sizeof(kmer_magic));
    stream->write(reinterpret_cast<char*>(&kmer_size), sizeof(kmer_size));
    stream->write(reinterpret_cast<char*>(&kmer_slots), sizeof(kmer_slots));
    stream->write(reinterpret_cast<char*>(&count_slots), sizeof(count_slots));
    stream->write(reinterpret_cast<char*>(&id), sizeof(id));
    stream->write(reinterpret_cast<char*>(&partition), sizeof(partition));
  }

  void deserialize(std::istream* stream) override
  {
    _deserialize(stream);
    stream->read(reinterpret_cast<char*>(&kmer_magic), sizeof(kmer_magic));
    stream->read(reinterpret_cast<char*>(&kmer_size), sizeof(kmer_size));
    stream->read(reinterpret_cast<char*>(&kmer_slots), sizeof(kmer_slots));
    stream->read(reinterpret_cast<char*>(&count_slots), sizeof(count_slots));
    stream->read(reinterpret_cast<char*>(&id), sizeof(id));
    stream->read(reinterpret_cast<char*>(&partition), sizeof(partition));
  }

  void sanity_check() override
  {
    _sanity_check();
    if (kmer_magic != MAGICS.at(KM_FILE::KMER))
      throw IOError("Invalid file format.");
  }

public:
  uint64_t kmer_magic {MAGICS.at(KM_FILE::KMER)};
  uint32_t kmer_size;
  uint32_t kmer_slots;
  uint32_t count_slots;
  uint32_t id;
  uint32_t partition;
};

template<size_t buf_size = 8192>
class KmerWriter : public IFile<KmerFileHeader, std::ostream, buf_size>
{
  using ocstream = lz4_stream::basic_ostream<buf_size>;
public:
  KmerWriter(const std::string& path,
             uint32_t kmer_size,
             uint32_t count_size,
             uint32_t id,
             uint32_t partition,
             bool lz4)
    : IFile<KmerFileHeader, std::ostream, buf_size>(path, std::ios::out | std::ios::binary)
  {
    this->m_header.compressed = lz4;
    this->m_header.kmer_size = kmer_size;
    this->m_header.kmer_slots = (kmer_size + 31) / 32;
    this->m_header.count_slots = count_size;
    this->m_header.id = id;
    this->m_header.partition = partition;

    this->m_header.serialize(this->m_first_layer.get());

    this->template set_second_layer<ocstream>(this->m_header.compressed);
  }

  template<size_t MAX_K, size_t MAX_C>
  void write(const Kmer<MAX_K>& kmer, const typename selectC<MAX_C>::type count)
  {
    this->m_second_layer->write(reinterpret_cast<const char*>(kmer.get_data64()),
                                this->m_header.kmer_slots*8);
    this->m_second_layer->write(reinterpret_cast<const char*>(&count), sizeof(count));
  }

  template<size_t MAX_C>
  void write_raw(const uint64_t* data, const typename selectC<MAX_C>::type count)
  {
    this->m_second_layer->write(reinterpret_cast<const char*>(data),
                                this->m_header.kmer_slots*8);
    this->m_second_layer->write(reinterpret_cast<const char*>(&count), sizeof(count));
  }
};

template<size_t buf_size>
using kw_t = std::shared_ptr<KmerWriter<buf_size>>;

template<size_t buf_size = 8192>
class KmerReader : public IFile<KmerFileHeader, std::istream, buf_size>
{
  using icstream = lz4_stream::basic_istream<buf_size>;
public:
  KmerReader(const std::string& path)
    : IFile<KmerFileHeader, std::istream, buf_size>(path, std::ios::in | std::ios::binary)
  {
    this->m_header.deserialize(this->m_first_layer.get());
    this->m_header.sanity_check();
    this->template set_second_layer<icstream>(this->m_header.compressed);
  }

  template<size_t MAX_K, size_t MAX_C>
  bool read(Kmer<MAX_K>& kmer, typename selectC<MAX_C>::type& count)
  {
    this->m_second_layer->read(reinterpret_cast<char*>(kmer.get_data64_unsafe()),
                               this->m_header.kmer_slots*8);
    this->m_second_layer->read(reinterpret_cast<char*>(&count), this->m_header.count_slots);

    if (!this->m_second_layer->gcount())
      return false;
    return true;
  }

  template<size_t MAX_K, size_t MAX_C>
  void write_as_text(std::ostream& stream)
  {
    Kmer<MAX_K> kmer; kmer.set_k(this->m_header.kmer_size);
    typename selectC<MAX_C>::type count = 0;
    while (read<MAX_K, MAX_C>(kmer, count))
    {
      stream << kmer.to_string() << " " << std::to_string(count) << "\n";
    }
  }

  template<size_t MAX_K, size_t MAX_C>
  void write_kmers(std::ostream& stream)
  {
    Kmer<MAX_K> kmer; kmer.set_k(this->m_header.kmer_size);
    typename selectC<MAX_C>::type count = 0;
    while (read<MAX_K, MAX_C>(kmer, count))
    {
      stream << kmer.to_string() << '\n';
    }
  }
};

template<size_t buf_size>
using kr_t = std::shared_ptr<KmerReader<buf_size>>;


template<size_t MAX_K, size_t MAX_C>
class KmerFileMerger
{
  using count_type = typename selectC<MAX_C>::type;
public:
  struct element
  {
    Kmer<MAX_K> value;
    count_type count {0};
    bool is_set {false};
  };

public:
  KmerFileMerger(const std::vector<std::string>& paths, uint32_t kmer_size)
    : m_paths(paths), m_kmer_size(kmer_size)
  {
    init_stream();
    init_state();
  }

  const Kmer<MAX_K>& current() const
  {
    return m_current;
  }

  count_type count() const
  {
    return m_counts;
  }

  void init_stream()
  {
    for (auto& path: m_paths)
      m_input_streams.push_back(std::make_shared<KmerReader<8192>>(path));
    m_size = m_paths.size();
  }

  void init_state()
  {
    for (size_t i=0; i<m_size; i++)
    {
      m_elements.push_back(element{});
      m_elements[i].value.set_k(m_kmer_size);

      if (read_next(i))
        m_elements[i].is_set = true;

      if ( (!m_current_set || m_elements[i].value < m_current) && m_elements[i].is_set )
      {
        m_current = m_next = m_elements[i].value;
        m_current_set = true;
      }
    }
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
    KmerWriter<8192> kw(path, m_kmer_size, requiredC<MAX_C>::value/8, 0, -1, compressed);
    while (next())
    {
      kw.template write<MAX_K, MAX_C>(m_current, m_counts);
    }
  }

  void write_as_text(std::ostream& out)
  {
    while (next())
    {
      out << m_current.to_string() << " " << std::to_string(m_counts) << "\n";
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

  std::vector<kr_t<8192>> m_input_streams;
  std::vector<element> m_elements;

  uint32_t m_size;
  uint32_t m_kmer_size;

  Kmer<MAX_K> m_next;
  Kmer<MAX_K> m_current;
  bool m_next_set {false};
  bool m_current_set {false};
  count_type m_counts;

  bool m_finish {false};
};

template<size_t MAX_K, size_t MAX_C>
class KmerFileAggregator
{
public:
  KmerFileAggregator(const std::vector<std::string>& paths, uint32_t kmer_size)
    : m_paths(paths), m_kmer_size(kmer_size)
  {

  }

  void write_as_bin(const std::string& path, bool compressed)
  {
    KmerWriter<8192> kw(path, m_kmer_size, requiredC<MAX_C>::value/8, 0, -1, compressed);
    for (auto& p : m_paths)
    {
      KmerReader<8192> kr(p);
      Kmer<MAX_K> k; k.set_k(m_kmer_size);
      typename selectC<MAX_C>::type count;
      while (kr.template read<MAX_K, MAX_C>(k, count))
      {
        kw.template write<MAX_K, MAX_C>(k, count);
      }
    }
  }

  void write_as_text(std::ostream& out)
  {
    for (auto& p : m_paths)
    {
      KmerReader<8192> kr(p);
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
      KmerReader<8192> kr(p);
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

};

