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

namespace km {

class SuperkFileHeader : public KmHeader
{
public:
  SuperkFileHeader() {};

  void serialize(std::ostream* stream) override
  {
    _serialize(stream);
    stream->write(reinterpret_cast<char*>(&superk_magic), sizeof(superk_magic));
    stream->write(reinterpret_cast<char*>(&partition), sizeof(partition));
  }

  void deserialize(std::istream* stream) override
  {
    _deserialize(stream);
    stream->read(reinterpret_cast<char*>(&superk_magic), sizeof(superk_magic));
    stream->read(reinterpret_cast<char*>(&partition), sizeof(partition));
  }

  void sanity_check() override
  {
    _sanity_check();
    if (superk_magic!= MAGICS.at(KM_FILE::SUPERK))
      throw IOError("Invalid file format.");
  }

public:
  uint64_t superk_magic {MAGICS.at(KM_FILE::SUPERK)};
  uint32_t partition;
};

template<size_t buf_size = 8192>
class SuperkWriter : public IFile<SuperkFileHeader, std::ostream, buf_size>
{
  using ocstream = lz4_stream::basic_ostream<buf_size>;
public:
  SuperkWriter(const std::string& path,
             uint32_t partition,
             bool lz4)
    : IFile<SuperkFileHeader, std::ostream, buf_size>(path, std::ios::out | std::ios::binary)
  {
    this->m_header.compressed = lz4;
    this->m_header.partition = partition;

    this->m_header.serialize(this->m_first_layer.get());

    this->template set_second_layer<ocstream>(this->m_header.compressed);
  }

  void write_size(unsigned int size)
  {
    this->m_second_layer->write(reinterpret_cast<char*>(&size), sizeof(size));
  }

  void write_block(unsigned char* block, size_t size)
  {
    this->m_second_layer->write(reinterpret_cast<char*>(block), size);
  }
};

template<size_t buf_size = 8192>
using skw_t = std::shared_ptr<SuperkWriter<buf_size>>;

template<size_t buf_size = 8192>
class SuperkReader : public IFile<SuperkFileHeader, std::istream, buf_size>
{
  using icstream = lz4_stream::basic_istream<buf_size>;
public:
  SuperkReader(const std::string& path)
    : IFile<SuperkFileHeader, std::istream, buf_size>(path, std::ios::in | std::ios::binary)
  {
    this->m_header.deserialize(this->m_first_layer.get());
    this->m_header.sanity_check();

    this->template set_second_layer<icstream>(this->m_header.compressed);
  }

  int read_size(unsigned int* size)
  {
    this->m_second_layer->read(reinterpret_cast<char*>(size), sizeof(*size));

    if (this->m_second_layer->gcount())
      return this->m_second_layer->gcount();
    return 0;
  }

  bool read_block(unsigned char* block, size_t size)
  {
    this->m_second_layer->read(reinterpret_cast<char*>(block), size);

    if (!this->m_second_layer->gcount())
      return false;
    return true;
  }
};

template<size_t buf_size = 8192>
using skr_t = std::shared_ptr<SuperkReader<buf_size>>;

};