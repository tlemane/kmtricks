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
#include <kmtricks/utils.hpp>
#include <vp4.h>

namespace km {

class HashFileHeader : public KmHeader
{
public:
  HashFileHeader() {};

  void serialize(std::ostream* stream) override
  {
    _serialize(stream);
    stream->write(reinterpret_cast<char*>(&hash_magic), sizeof(hash_magic));
    stream->write(reinterpret_cast<char*>(&count_slots), sizeof(count_slots));
    stream->write(reinterpret_cast<char*>(&id), sizeof(id));
    stream->write(reinterpret_cast<char*>(&partition), sizeof(partition));
  }

  void deserialize(std::istream* stream) override
  {
    _deserialize(stream);
    stream->read(reinterpret_cast<char*>(&hash_magic), sizeof(hash_magic));
    stream->read(reinterpret_cast<char*>(&count_slots), sizeof(count_slots));
    stream->read(reinterpret_cast<char*>(&id), sizeof(id));
    stream->read(reinterpret_cast<char*>(&partition), sizeof(partition));
  }

  void sanity_check() override
  {
    _sanity_check();
    if (hash_magic != MAGICS.at(KM_FILE::HASH))
      throw IOError("Invalid file format.");
  }

public:
  uint64_t hash_magic {MAGICS.at(KM_FILE::HASH)};
  uint32_t count_slots;
  uint32_t id;
  uint32_t partition;
};

template<size_t MAX_C, size_t buf_size = 32768>
class HashWriter : public IFile<HashFileHeader, std::ostream, buf_size>
{
  using ocstream = lz4_stream::basic_ostream<buf_size>;
  using count_type = typename selectC<MAX_C>::type;
public:
  HashWriter(const std::string& path,
             uint32_t count_size,
             uint32_t id,
             uint32_t partition,
             bool compress)
    : IFile<HashFileHeader, std::ostream, buf_size>(path, std::ios::out | std::ios::binary)
  {
    this->m_header.compressed = compress;
    this->m_header.count_slots = count_size;
    this->m_header.id = id;
    this->m_header.partition = partition;

    this->m_header.serialize(this->m_first_layer.get());
    this->template set_second_layer<ocstream>(false);
  }

  ~HashWriter()
  {
    if (m_index != 0)
      flush();
  }

  void write(uint64_t hash, count_type count)
  {
    if (m_index == m_capacity)
      flush();
    m_src[m_index] = hash;
    m_src_c[m_index] = count;
    m_index++;
  }

  void flush()
  {
    if (!m_index)
      return;

    if (this->m_header.compressed)
    {
      size_t n = m_index;
      size_t hash_bytes = p4nd1enc64(m_src.data(), n, m_dest.data());
      size_t count_bytes = 0;
      if constexpr(sizeof(count_type) == 1)
        count_bytes = p4nzenc8(m_src_c.data(), n, m_dest_c.data());
      else if constexpr(sizeof(count_type) == 2)
        count_bytes = p4nzenc16(m_src_c.data(), n, m_dest_c.data());
      else
        count_bytes = p4nzenc32(m_src_c.data(), n, m_dest_c.data());

      this->m_second_layer->write(reinterpret_cast<char*>(&n), sizeof(n));
      this->m_second_layer->write(reinterpret_cast<char*>(&hash_bytes), sizeof(hash_bytes));
      this->m_second_layer->write(reinterpret_cast<char*>(m_dest.data()), hash_bytes);
      this->m_second_layer->write(reinterpret_cast<char*>(&count_bytes), sizeof(count_bytes));
      this->m_second_layer->write(reinterpret_cast<char*>(m_dest_c.data()), count_bytes);

    }
    else
    {
      this->m_second_layer->write(reinterpret_cast<char*>(&m_index), sizeof(m_index));
      this->m_second_layer->write(reinterpret_cast<char*>(m_src.data()), m_index * sizeof(uint64_t));
      this->m_second_layer->write(reinterpret_cast<char*>(m_src_c.data()), m_index * sizeof(count_type));
    }
    m_index = 0;
  }

private:
  std::array<uint64_t, buf_size / sizeof(uint64_t)> m_src;
  std::array<unsigned char, buf_size> m_dest;
  std::array<count_type, buf_size / sizeof(count_type)> m_src_c;
  std::array<unsigned char, buf_size> m_dest_c;
  size_t m_index {0};
  size_t m_in_buffer {0};
  size_t m_capacity {m_src.size()};
};

template<size_t MAX_C, size_t buf_size = 32768>
using hw_t = std::shared_ptr<HashWriter<MAX_C, buf_size>>;

template<size_t MAX_C, size_t buf_size = 32768>
class HashReader : public IFile<HashFileHeader, std::istream, buf_size>
{
  using icstream = lz4_stream::basic_istream<buf_size>;
  using count_type = typename selectC<MAX_C>::type;

public:
  HashReader(const std::string& path)
    : IFile<HashFileHeader, std::istream, buf_size>(path, std::ios::in | std::ios::binary)
  {
    this->m_header.deserialize(this->m_first_layer.get());
    this->m_header.sanity_check();

    this->template set_second_layer<icstream>(false);
  }

  bool load()
  {
    this->m_second_layer->read(reinterpret_cast<char*>(&m_in_buffer), sizeof(m_in_buffer));
    if (!this->m_second_layer->gcount())
      return false;

    if (this->m_header.compressed)
    {
      size_t hash_bytes = 0;
      size_t count_bytes = 0;

      this->m_second_layer->read(reinterpret_cast<char*>(&hash_bytes), sizeof(hash_bytes));
      this->m_second_layer->read(reinterpret_cast<char*>(m_src.data()), hash_bytes);
      this->m_second_layer->read(reinterpret_cast<char*>(&count_bytes), sizeof(count_bytes));
      this->m_second_layer->read(reinterpret_cast<char*>(m_src_c.data()), count_bytes);

      p4nd1dec64(m_src.data(), m_in_buffer, m_dest.data());

      if constexpr(sizeof(count_type) == 1)
        p4nzdec8(m_src_c.data(), m_in_buffer, m_dest_c.data());
      else if constexpr(sizeof(count_type) == 2)
        p4nzdec16(m_src_c.data(), m_in_buffer, m_dest_c.data());
      else
        p4nzdec32(m_src_c.data(), m_in_buffer, m_dest_c.data());
    }
    else
    {
      this->m_second_layer->read(reinterpret_cast<char*>(m_dest.data()), m_in_buffer * sizeof(uint64_t));
      this->m_second_layer->read(reinterpret_cast<char*>(m_dest_c.data()), m_in_buffer * sizeof(count_type));
    }
    m_index = 0;
    return true;
  }

  bool read(uint64_t& hash, count_type& count)
  {
    if (m_in_buffer == 0)
      if (!load())
        return false;

    hash = m_dest[m_index];
    count = m_dest_c[m_index];

    m_in_buffer--;
    m_index++;

    return true;
  }

  void write_as_text(std::ostream& stream)
  {
    uint64_t hash = 0;
    count_type count = 0;
    while (read(hash, count))
    {
      stream << std::to_string(hash) << " " << std::to_string(count) << "\n";
    }
  }

private:
  std::array<uint64_t, buf_size / sizeof(uint64_t)> m_dest;
  std::array<unsigned char, buf_size> m_src;
  std::array<count_type, buf_size / sizeof(count_type)> m_dest_c;
  std::array<unsigned char, buf_size> m_src_c;
  size_t m_index {0};
  size_t m_in_buffer {0};
  size_t m_capacity {m_dest.size()};
};

template<size_t MAX_C, size_t buf_size = 32768>
using hr_t = std::shared_ptr<HashReader<MAX_C, buf_size>>;

template<size_t MAX_C>
class HashFileAggregator
{
  using count_type = typename selectC<MAX_C>::type;
public:
  HashFileAggregator(const std::vector<std::string>& paths)
    : m_paths(paths)
  {

  }

  void write_as_bin(const std::string& path, bool compressed)
  {
    HashWriter<MAX_C, 32768> kw(path, requiredC<MAX_C>::value/8, 0, -1, compressed);
    for (auto& p : m_paths)
    {
      HashReader<MAX_C, 32768> kr(p);
      uint64_t hash;
      count_type count;
      while (kr.read(hash, count))
      {
        kw.write(hash, count);
      }
    }
  }

  void write_as_text(std::ostream& out)
  {
    for (auto& p : m_paths)
    {
      HashReader<MAX_C, 32768> kr(p);
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