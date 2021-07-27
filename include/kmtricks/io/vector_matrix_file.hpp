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
#include <kmtricks/bitmatrix.hpp>

namespace km {

class VectorMatrixFileHeader : public KmHeader
{
public:
  VectorMatrixFileHeader() {};

  void serialize(std::ostream* stream)
  {
    _serialize(stream);
    stream->write(reinterpret_cast<char*>(&matrix_magic), sizeof(matrix_magic));
    stream->write(reinterpret_cast<char*>(&bits), sizeof(bits));
    stream->write(reinterpret_cast<char*>(&first), sizeof(first));
    stream->write(reinterpret_cast<char*>(&window), sizeof(window));
    stream->write(reinterpret_cast<char*>(&id), sizeof(id));
    stream->write(reinterpret_cast<char*>(&partition), sizeof(partition));
  }

  void deserialize(std::istream* stream)
  {
    _deserialize(stream);
    stream->read(reinterpret_cast<char*>(&matrix_magic), sizeof(matrix_magic));
    stream->read(reinterpret_cast<char*>(&bits), sizeof(bits));
    stream->read(reinterpret_cast<char*>(&first), sizeof(first));
    stream->read(reinterpret_cast<char*>(&window), sizeof(window));
    stream->read(reinterpret_cast<char*>(&id), sizeof(id));
    stream->read(reinterpret_cast<char*>(&partition), sizeof(partition));
  }

  void sanity_check()
  {
    _sanity_check();
    if (matrix_magic != MAGICS.at(KM_FILE::BITMATRIX))
      throw IOError("Invalid file format.");
  }

public:
  uint64_t matrix_magic {MAGICS.at(KM_FILE::BITMATRIX)};
  uint32_t bits;
  uint32_t id;
  uint32_t partition;
  uint64_t first;
  uint64_t window;
};

template<size_t buf_size = 8192>
class VectorMatrixWriter : public IFile<VectorMatrixFileHeader, std::ostream, buf_size>
{
  using ocstream = lz4_stream::basic_ostream<buf_size>;
public:
  VectorMatrixWriter(const std::string& path,
                   uint32_t bits,
                   uint32_t id,
                   uint32_t partition,
                   uint64_t first,
                   uint64_t window,
                   bool lz4)
    : IFile<VectorMatrixFileHeader, std::ostream, buf_size>(path, std::ios::out | std::ios::binary)
  {
    this->m_header.compressed = lz4;
    this->m_header.bits = bits;
    this->m_header.first = first;
    this->m_header.window = window;
    this->m_header.id = id;
    this->m_header.partition = partition;

    this->m_header.serialize(this->m_first_layer.get());

    this->template set_second_layer<ocstream>(this->m_header.compressed);
  }

  void write(std::vector<uint8_t>& bits)
  {
    this->m_second_layer->write(reinterpret_cast<char*>(bits.data()), bits.size()*sizeof(uint8_t));
  }

  void dump(BitMatrix& bit_matrix)
  {
    this->m_second_layer->write(reinterpret_cast<char*>(bit_matrix.matrix),
                                bit_matrix.get_size_in_byte());
  }
};

template<size_t buf_size = 8192>
class VectorMatrixReader : public IFile<VectorMatrixFileHeader, std::istream, buf_size>
{
  using icstream = lz4_stream::basic_istream<buf_size>;
public:
  VectorMatrixReader(const std::string& path)
    : IFile<VectorMatrixFileHeader, std::istream, buf_size>(path, std::ios::in | std::ios::binary)
  {
    this->m_header.deserialize(this->m_first_layer.get());
    this->m_header.sanity_check();

    if (this->m_header.compressed)
      this->template set_second_layer<icstream>(this->m_header.compressed);
  }

  bool read(std::vector<uint8_t>& bits)
  {
    if (this->m_header.compressed)
    {
      this->m_second_layer->read(reinterpret_cast<char*>(bits.data()), bits.size()*sizeof(uint8_t));
      if (!this->m_second_layer->gcount())
        return false;
    }
    else
    {
      this->m_first_layer->read(reinterpret_cast<char*>(bits.data()), bits.size()*sizeof(uint8_t));
      if (!this->m_first_layer->gcount())
        return false;
    }
    return true;
  }

  bool read(char* bits, size_t size)
  {
    if (this->m_header.compressed)
    {
      this->m_second_layer->read(bits, size);
      if (!this->m_second_layer->gcount())
        return false;
    }
    else
    {
      this->m_first_layer->read(bits, size);
      if (!this->m_first_layer->gcount())
        return false;
    }
    return true;
  }

  void load(BitMatrix& bit_matrix)
  {
    if (this->m_header.compressed)
    {
      this->m_second_layer->read(reinterpret_cast<char*>(bit_matrix.matrix),
                                 bit_matrix.get_size_in_byte());
    }
    else
    {
      this->m_first_layer->read(reinterpret_cast<char*>(bit_matrix.matrix),
                                bit_matrix.get_size_in_byte());
    }
  }

  void seekg(uint32_t partition)
  {
    if (this->m_header.compressed)
      throw IOError("VectorMatrixReader::seekg() only available on uncompressed stream.");
    this->m_first_layer->seekg(49 + (partition * (this->m_header.window / 8)));
  }
};

template<size_t buf_size = 8192>
using vmr_t = std::shared_ptr<VectorMatrixReader<buf_size>>;

};