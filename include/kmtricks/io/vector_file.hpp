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

class BitVectorFileHeader : public KmHeader
{
public:
  BitVectorFileHeader() {};

  void serialize(std::ostream* stream)
  {
    _serialize(stream);
    stream->write(reinterpret_cast<char*>(&bit_vector_magic), sizeof(bit_vector_magic));
    stream->write(reinterpret_cast<char*>(&bits), sizeof(bits));
    stream->write(reinterpret_cast<char*>(&id), sizeof(id));
    stream->write(reinterpret_cast<char*>(&partition), sizeof(partition));
  }

  void deserialize(std::istream* stream)
  {
    _deserialize(stream);
    stream->read(reinterpret_cast<char*>(&bit_vector_magic), sizeof(bit_vector_magic));
    stream->read(reinterpret_cast<char*>(&bits), sizeof(bits));
    stream->read(reinterpret_cast<char*>(&id), sizeof(id));
    stream->read(reinterpret_cast<char*>(&partition), sizeof(partition));
  }

  void sanity_check()
  {
    _sanity_check();
    if (bit_vector_magic != MAGICS.at(KM_FILE::VECTOR))
      throw IOError("Invalid file format.");
  }

public:
  uint64_t bit_vector_magic {MAGICS.at(KM_FILE::VECTOR)};
  uint64_t bits;
  uint32_t id;
  uint32_t partition;
};

template<size_t buf_size = 8192>
class BitVectorWriter : public IFile<BitVectorFileHeader, std::ostream, buf_size>
{
  using ocstream = lz4_stream::basic_ostream<buf_size>;
public:
  BitVectorWriter(const std::string& path, uint64_t bits, uint32_t id, uint32_t partition, bool lz4)
    : IFile<BitVectorFileHeader, std::ostream, buf_size>(path, std::ios::out | std::ios::binary)
  {
    this->m_header.compressed = lz4;
    this->m_header.bits = bits;
    this->m_header.id = id;
    this->m_header.partition = partition;

    this->m_header.serialize(this->m_first_layer.get());

    this->template set_second_layer<ocstream>(this->m_header.compressed);
  }

  void write(std::vector<uint8_t>& bits)
  {
    this->m_second_layer->write(reinterpret_cast<char*>(bits.data()), bits.size()*sizeof(uint8_t));
  }
};

template<size_t buf_size = 8192>
using bvw_t = std::shared_ptr<BitVectorWriter<buf_size>>;

template<size_t buf_size = 8192>
class BitVectorReader : public IFile<BitVectorFileHeader, std::istream, buf_size>
{
  using icstream = lz4_stream::basic_istream<buf_size>;
public:
  BitVectorReader(const std::string& path)
    : IFile<BitVectorFileHeader, std::istream, buf_size>(path, std::ios::in | std::ios::binary)
  {
    this->m_header.deserialize(this->m_first_layer.get());
    this->m_header.sanity_check();

    this->template set_second_layer<icstream>(this->m_header.compressed);
  }

  void read(std::vector<uint8_t>& bits)
  {
    this->m_second_layer->read(reinterpret_cast<char*>(bits.data()), bits.size()*sizeof(uint8_t));
  }

  void read(char* bits, size_t size)
  {
    this->m_second_layer->read(bits, size);
  }
};

};