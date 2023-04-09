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
#include <string>
#include <fstream>
#include <array>
#include <memory>
#include <map>

#include <kmtricks/io/lz4_stream.hpp>
#include <kmtricks/exceptions.hpp>
#include <kmtricks/utils.hpp>

#define KM_IO_VERSION 0x0

namespace km {

enum class KM_FILE
{
  BASE,
  KMER,
  HASH,
  MATRIX,
  MATRIX_HASH,
  PAMATRIX,
  PAMATRIX_HASH,
  VECTOR,
  BITMATRIX,
  KFF,
  HIST,
  SUPERK
};

const std::map<KM_FILE, uint64_t> MAGICS = {
  {KM_FILE::BASE, 0x736b636972746d6b},
  {KM_FILE::KMER, 0x72656d6b},
  {KM_FILE::HASH, 0x68736168},
  {KM_FILE::MATRIX, 0x6b5f78697274616d},
  {KM_FILE::PAMATRIX, 0x6b5f74616d6170},
  {KM_FILE::VECTOR, 0x726f74636576},
  {KM_FILE::BITMATRIX, 0x74616d746962},
  {KM_FILE::HIST, 0x747369686b},
  {KM_FILE::SUPERK, 0x6b7265707573},
  {KM_FILE::MATRIX_HASH, 0x685f78697274616d},
  {KM_FILE::PAMATRIX_HASH, 0x685f74616d6170}
};

inline KM_FILE get_km_file_type(const std::string& path)
{
  uint64_t km_base;
  uint64_t km_file;
  std::ifstream in(path, std::ios::in | std::ios::binary); check_fstream_good(path, in);
  in.read(reinterpret_cast<char*>(&km_base), sizeof(km_base));
  if (km_base != MAGICS.at(KM_FILE::BASE))
    throw IOError("Not a kmtricks file.");
  in.ignore(5);
  in.read(reinterpret_cast<char*>(&km_file), sizeof(km_file));

  if (km_file == MAGICS.at(KM_FILE::KMER))
    return KM_FILE::KMER;
  else if (km_file == MAGICS.at(KM_FILE::HASH))
    return KM_FILE::HASH;
  else if (km_file == MAGICS.at(KM_FILE::MATRIX))
    return KM_FILE::MATRIX;
  else if (km_file == MAGICS.at(KM_FILE::MATRIX_HASH))
    return KM_FILE::MATRIX_HASH;
  else if (km_file == MAGICS.at(KM_FILE::PAMATRIX))
    return KM_FILE::PAMATRIX;
  else if (km_file == MAGICS.at(KM_FILE::PAMATRIX_HASH))
    return KM_FILE::PAMATRIX_HASH;
  else if (km_file == MAGICS.at(KM_FILE::VECTOR))
    return KM_FILE::VECTOR;
  else if (km_file == MAGICS.at(KM_FILE::BITMATRIX))
    return KM_FILE::BITMATRIX;
  else if (km_file == MAGICS.at(KM_FILE::HIST))
    return KM_FILE::HIST;
  else if (km_file == MAGICS.at(KM_FILE::SUPERK))
    return KM_FILE::SUPERK;
  else
    throw IOError("Not a kmtricks file.");
}

inline std::string km_file_to_str(KM_FILE f)
{
  if (f == KM_FILE::KMER)
    return "kmer";
  else if (f == KM_FILE::HASH)
    return "hash";
  else if (f == KM_FILE::MATRIX)
    return "count matrix";
  else if (f == KM_FILE::MATRIX_HASH)
    return "hash matrix";
  else if (f == KM_FILE::PAMATRIX)
    return "pa matrix";
  else if (f == KM_FILE::PAMATRIX_HASH)
    return "hash pa matrix";
  else if (f == KM_FILE::VECTOR)
    return "bit vector";
  else if (f == KM_FILE::BITMATRIX)
    return "bit matrix";
  else if (f == KM_FILE::HIST)
    return "histogram";
  else if (f == KM_FILE::SUPERK)
    return "super-k-mer";
  else
    return "base";
}

class KmHeader
{
public:
  KmHeader() {}
  virtual void sanity_check() = 0;
  virtual void serialize(std::ostream*) = 0;
  virtual void deserialize(std::istream*) = 0;

protected:
  void _serialize(std::ostream* stream)
  {
    stream->write(reinterpret_cast<char*>(&km_magic), sizeof(km_magic));
    stream->write(reinterpret_cast<char*>(&km_version), sizeof(km_version));
    stream->write(reinterpret_cast<char*>(&compressed), sizeof(compressed));
  }

  void _deserialize(std::istream* stream)
  {
    stream->read(reinterpret_cast<char*>(&km_magic), sizeof(km_magic));
    stream->read(reinterpret_cast<char*>(&km_version), sizeof(km_version));
    stream->read(reinterpret_cast<char*>(&compressed), sizeof(compressed));
  }

  void _sanity_check()
  {
    if (km_magic != MAGICS.at(KM_FILE::BASE))
      throw IOError("Invalid file format.");
  }

public:
  uint64_t km_magic {MAGICS.at(KM_FILE::BASE)};
  uint32_t km_version {KM_IO_VERSION};
  bool compressed;
};

template<typename header_t,
         typename stream,
         size_t   buf_size>
class IFile
{
  using stream_t = std::unique_ptr<stream>;
  using buffer_t = std::array<char, buf_size>;

public:
  IFile () : m_first_layer(new std::fstream{}) {}

  IFile (const std::string& path, std::ios_base::openmode mode)
    : m_first_layer(new std::fstream{path, mode}), m_path(path)
  {
    if (!this->m_first_layer->good())
      throw std::runtime_error("Unable to open " + path);
  }
  IFile (IFile const &) = delete;
  IFile& operator= (IFile const &) = delete;
  IFile (IFile&&) = delete;
  IFile& operator= (IFile&&) = delete;

  virtual ~IFile() = default;

  const header_t& infos() const
  {
    return m_header;
  }

  template<typename compression_stream_t>
  stream_t add_compression_layer(bool compressed)
  {
    if (compressed)
      return std::unique_ptr<compression_stream_t>(
        new compression_stream_t(*this->m_first_layer.get())
      );
    else
      return std::move(this->m_first_layer);
  }

  stream* get_stream()
  {
    return m_second_layer->get();
  }

  void flush()
  {
  }

  void close()
  {
  }

protected:

  template<typename compression_stream_t>
  void set_second_layer(bool compress)
  {
    this->m_second_layer = this->template add_compression_layer<compression_stream_t>(compress);

#ifndef __APPLE__
    // Dirty fix, I don't know why but this line cause a segfault on Apple clang > 10.1
    this->m_second_layer->rdbuf()->pubsetbuf(m_buf.data(), m_buf.size());
#endif
  }

protected:
  stream_t    m_first_layer  {nullptr}; // fstream layer
  stream_t    m_second_layer {nullptr}; // compression layer, must inherit from basic_xstream<char>
  buffer_t    m_buf;
  header_t    m_header;
  std::string m_path;
};

};
