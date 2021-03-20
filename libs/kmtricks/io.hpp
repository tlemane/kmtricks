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
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <array>
#include <memory>

#include "sequences.hpp"
#include "lz4_stream.hpp"
#include "utilities.hpp"
#include "bitmatrix.hpp"


//! \defgroup IO

//! \namespace km
//! The kmtricks library namespace
namespace km 
{

typedef std::basic_ostream<char> os;
typedef std::basic_istream<char> is;

using IN = is;
using OUT = os;

static const uint64_t magic1 = 0x41EFD2;
static const uint64_t magic2 = 0x03627B0E;

//! \ingroup IO
//! \brief enum class for matrix format
enum class matrix_format_t : uint32_t {
  ASCII, // <kmer> <count> (as text)
  BIN,   // <kmer> <count> (as uintXX_t)
  PA,    // <kmer> <bit-vector>
  BF,    // <bit-vector> (hash is index)
  BIT    // Only to store and load km::BitMatrix
};

typedef matrix_format_t matrix_t;

//! \ingroup IO
//! \brief map matrix_format_t to string
const static std::map<matrix_format_t, std::string> fmt_to_string = {
  {matrix_format_t::ASCII, "ASCII"},
  {matrix_format_t::BIN, "BIN"},
  {matrix_format_t::PA, "PA"},
  {matrix_format_t::BF, "BF"},
  {matrix_format_t::BIT, "BIT"}
};

//! \ingroup IO
//! \brief map string to matrix_format_t
const static std::map<std::string, matrix_format_t> string_to_fmt = {
  {"ASCII", matrix_format_t::ASCII},
  {"BIN", matrix_format_t::BIN},
  {"PA", matrix_format_t::PA},
  {"BF", matrix_format_t::BF},
  {"BIT", matrix_format_t::BIT}
};

struct kmer_file_header {
  uint64_t first_magic;
  uint32_t ktsize;        // bytes per kmer
  uint32_t ctsize;        // bytes per count
  uint32_t file_id;       // in kmtricks pipeline is the fof index
  uint32_t partition_id;  // partition id in [0..P-1]
  uint32_t kmer_size;     // size of a kmer
  uint32_t is_compressed; // 1 if compressed
  uint32_t is_hashes;     // 1 if hashes
  uint64_t second_magic;
};

struct count_matrix_file_header {
  uint64_t first_magic;
  matrix_t matrix_type;    // matrix_t::ASCII | matrix_t::BIN
  uint32_t ktsize;         // bytes per kmer
  uint32_t ctsize;         // bytes per count
  int32_t  partition_id;   // partition id
  uint32_t kmer_size;      // size of kmer
  uint32_t nb_counts;      // number of counts per kmer == nb samples
  uint32_t is_hashes;      // 1 if hashes
  uint32_t is_compressed;  // 1 if compressed
  uint64_t second_magic;
};

struct pa_matrix_file_header {
  uint64_t first_magic;
  matrix_t matrix_type;    // matrix_t::PA
  uint32_t ktsize;         // bytes per k-mers
  int32_t  partition_id;   // partition id
  uint32_t kmer_size;      // size of k-mer
  uint32_t bits_in_use;    // number of bits actually used in pa vector, because of padding
  uint32_t size_in_bytes;  // the size of one row (preence/absence bit-vector) in byte
  uint32_t is_hashes;      // 1 if hashes
  uint32_t is_compressed;  // 1 if compressed
  uint64_t second_magic;
}; 

struct bitvector_file_header {
  uint64_t first_magic;
  uint32_t file_id;         // in kmtricks pipeline is the fof index
  uint32_t partition_id;    // partition id in [0..P-1]
  uint32_t partition_size;  // partition size in bytes
  uint32_t is_compressed;   // 1 if compressed
  uint64_t nb_bits;         // bit vector length
  uint64_t second_magic;
};

struct bitmatrix_file_header {
  uint64_t first_magic;
  matrix_t matrix_type;   // matrix_t::BF, matrix_t::BIT
  int32_t  partition_id;  // partition id
  uint64_t nb_rows;       // number of rows
  uint64_t nb_rows_use;   // because padding
  uint64_t nb_cols;       // number of columns
  uint64_t nb_cols_use;   // because padding
  uint32_t size_in_bytes; // the size of one row in byte
  uint32_t is_compressed; // 1 if compressed
  uint64_t second_magic;
}; 

struct hist_file_header {
  uint64_t first_magic;
  int32_t  id;
  uint32_t kmer_size;
  uint64_t lower;
  uint64_t upper;
  uint64_t uniq;
  uint64_t total;
  uint64_t second_magic;
};

typedef struct bitmatrix_file_header bmheader_t;
typedef struct bitvector_file_header bvheader_t;
typedef struct pa_matrix_file_header paheader_t;
typedef struct count_matrix_file_header cmheader_t;
typedef struct kmer_file_header kheader_t;
typedef struct hist_file_header hheader_t;

const static kheader_t  kheader_d  = {magic1, 0, 0, 0, 0, 0, 0, 0, magic2};
const static cmheader_t cmheader_d = {magic1, matrix_t::ASCII, 0, 0, 0, 0, 0, 0, 0, magic2};
const static paheader_t paheader_d = {magic1, matrix_t::PA, 0, 0, 0, 0, 0, 0, 0, magic2};
const static bvheader_t bvheader_d = {magic1, 0, 0, 0, 0, 0, magic2};
const static bmheader_t bmheader_d = {magic1, matrix_t::BF, 0, 0, 0, 0, 0, 0, 0, magic2};
const static hheader_t  hheader_d  = {magic1, 0, 0, 0, 0, 0, 0, magic2};

//! \ingroup IO
//! \class IFile
//! \brief IFile interface
template<typename header_t,
         typename stream,
         size_t   buf_size>
class IFile
{
  using stream_t = std::unique_ptr<stream>;
  using buffer_t = std::array<char, buf_size>;

public:
  IFile () : first_layer(new std::fstream{}) {}
  
  IFile (const string& path, header_t default_header, std::ios_base::openmode mode)
    : first_layer(new std::fstream{path, mode}), header(default_header), path(path)
  {
    if (!this->first_layer->good())
      throw std::runtime_error("Unable to open " + path);
    memset(&header, 0, sizeof(header));
    header.first_magic = magic1;
    header.second_magic = magic2;
  }
  IFile (IFile const &) = delete;
  IFile& operator= (IFile const &) = delete;
  IFile (IFile&&) = delete;
  IFile& operator= (IFile&&) = delete;


  const header_t* infos() const
  {
    return &header;
  }

  //! \brief add a stream compressor
  template<typename compression_stream_t>
  stream_t add_compression_layer(bool compressed)
  {
    if (compressed)
      return std::unique_ptr<compression_stream_t>(
        new compression_stream_t(*first_layer.get())
      );
    else
      return std::move(this->first_layer);
  }

protected:

  template<typename compression_stream_t>
  void set_second_layer(bool compress)
  {
    this->second_layer = this->template add_compression_layer<compression_stream_t>(compress);

#ifndef __APPLE__
    // Dirty fix, I don't know why but this line cause a segfault on Apple clang > 10.1
    this->second_layer->rdbuf()->pubsetbuf(buf.data(), buf.size());
#endif
  }

protected:
  stream_t    first_layer  {nullptr}; // fstream layer
  stream_t    second_layer {nullptr}; // compression layer, must inherit from basic_xstream<char>
  buffer_t    buf;
  header_t    header;
  std::string path;
};

//! \ingroup IO
//! \class KmerFile
//! \brief TODO
template<typename stream,
         typename KmerType,
         typename CountType,
         size_t   buf_size = 4096>
         //typename = typename std::enable_if<std::is_integral<KmerType>::value &&
         //                                    std::is_integral<CountType>::value, void>::type>
class KmerFile : public IFile<kheader_t, stream, buf_size>
{

  using ocstream  = lz4_stream::basic_ostream<buf_size>;
  using icstream  = lz4_stream::basic_istream<buf_size>;

public:
  KmerFile () = delete;
  KmerFile (KmerFile const &) = delete;
  KmerFile& operator= (KmerFile const &) = delete;
  KmerFile (KmerFile&&) = delete;
  KmerFile& operator= (KmerFile&&) = delete;

  ~KmerFile() = default;

  //! \brief Constructor, read a kmer file (output of km_superk_to_kmer_counts module)
  //! \param path : path to the file
  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, is>{}, S>::type>
  KmerFile(const std::string& path) 
    : IFile<kheader_t, stream, buf_size>(
        path, kheader_d, std::ios::binary | std::ios::in)
  {
    this->first_layer->read(reinterpret_cast<char*>(&this->header), sizeof(this->header));
    if (this->header.first_magic != magic1 || this->header.second_magic != magic2)
      throw std::runtime_error("Unable to read " + this->path + ". Possibly due to bad format.");
    this->template set_second_layer<icstream>(this->header.is_compressed);
  }

  //! \brief Constructor, write a kmer file
  //! \param path : path to the file
  //! \param file_id : a file id, in kmtricks the fof index
  //! \param partition_id : partition id
  //! \param kmer_size : size of k-mers
  //! \param is_compressed : compress with lz4 (0/1)
  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, os>{}, S>::type>
  KmerFile(const std::string& path, uint32_t file_id,
           uint32_t partition_id, uint32_t kmer_size,
           uint32_t is_hashes, uint32_t is_compressed) 
    : IFile<kheader_t, stream, buf_size>(
        path, kheader_d, std::ios::binary | std::ios::out)
  {
    this->header.file_id       = file_id;
    this->header.partition_id  = partition_id;
    this->header.kmer_size     = kmer_size;
    this->header.is_compressed = is_compressed;
    this->header.is_hashes     = is_hashes;
    this->header.ktsize        = sizeof(KmerType);
    this->header.ctsize        = sizeof(CountType);
    
    this->first_layer->write(reinterpret_cast<char*>(&this->header), sizeof(this->header)); 
    this->first_layer->flush();
    this->template set_second_layer<ocstream>(this->header.is_compressed);
  }

  //! \brief write a kmer and its count
  //! \param kmer : reference to a Kmer<KmerType>
  //! \param count : reference to a Countype (integral type)
  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, os>{}, void>::type>
  void write(Kmer<KmerType>& kmer, CountType& count)
  {
    KmerType wtmp = kmer.value();
    this->second_layer->write(reinterpret_cast<char*>(&wtmp), sizeof(KmerType));
    this->second_layer->write(reinterpret_cast<char*>(&count), sizeof(CountType));
  }

  //! \brief write a kmer and its count
  //! \param kmer : reference to a KmerType (integral type)
  //! \param count : reference to a CountType (integral type)
  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, os>{}, void>::type>
  void write(KmerType& kmer, CountType& count)
  {
    this->second_layer->write(reinterpret_cast<char*>(&kmer), sizeof(KmerType));
    this->second_layer->write(reinterpret_cast<char*>(&count), sizeof(CountType));
  }

  //! \brief read a kmer and its count
  //! \param : reference to a Kmer<KmerType>
  //! \param : reference to a CountType (integral type)
  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, is>{}, void>::type>
  bool read(Kmer<KmerType>& kmer, CountType& count)
  {
    this->second_layer->read(reinterpret_cast<char*>(&tmp), sizeof(KmerType));
    if (this->second_layer->eof()) return false;
    this->second_layer->read(reinterpret_cast<char*>(&count), sizeof(CountType));
    kmer.set_kmer(tmp, this->header.kmer_size);
    return true;
  }

  //! \brief read a kmer and its count
  //! \param kmer : reference to a KmerType (integral type)
  //! \param count : reference to CountType (integral type)
  //! \return false if end of file, else true
  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, is>{}, void>::type>
  bool read(KmerType& kmer, CountType& count)
  {
    this->second_layer->read(reinterpret_cast<char*>(&kmer), sizeof(KmerType));
    if (this->second_layer->eof()) return false;
    this->second_layer->read(reinterpret_cast<char*>(&count), sizeof(CountType));
    return true;
  }
private:
  KmerType tmp {0};
};

//! \ingroup IO
//! \class CountMatrix
//! \brief TODO
template<typename stream,
         typename KmerType,
         typename CountType,
         matrix_format_t MatrixType = matrix_format_t::ASCII,
         size_t buf_size = 4096,
         typename = typename std::enable_if<(std::is_integral<KmerType>::value ||
                                            std::is_same<KmerType, __uint128_t>{}) &&
                                            std::is_integral<CountType>::value &&
                                            (MatrixType == matrix_format_t::ASCII ||
                                            MatrixType == matrix_format_t::BIN), void>::type>
class CountMatrixFile : public IFile<cmheader_t, stream, buf_size>
{

  using ocstream = lz4_stream::basic_ostream<buf_size>;
  using icstream = lz4_stream::basic_istream<buf_size>;

private:
  void write_ascii_header()
  {
    *this->first_layer.get() << "FM " << std::to_string(this->header.first_magic) << "\n";
    *this->first_layer.get() << "MT " << fmt_to_string.at(this->header.matrix_type) << "\n";
    *this->first_layer.get() << "PT " << std::to_string(this->header.partition_id) << "\n";
    *this->first_layer.get() << "KS " << std::to_string(this->header.kmer_size) << "\n";
    *this->first_layer.get() << "NC " << std::to_string(this->header.nb_counts) << "\n";
    *this->first_layer.get() << "IH " << std::to_string(this->header.is_hashes) << "\n";
    *this->first_layer.get() << "IC " << std::to_string(this->header.is_compressed) << "\n";
    *this->first_layer.get() << "SM " << std::to_string(this->header.second_magic) << "\n";
  }
  void read_ascii_header()
  {
    string line;
    getline(*this->first_layer.get(), line); this->header.first_magic = stoi(split(line, ' ')[1]);
    getline(*this->first_layer.get(), line); this->header.matrix_type = string_to_fmt.at(
      split(line, ' ')[1]);
    getline(*this->first_layer.get(), line); this->header.partition_id = stoi(split(line, ' ')[1]);
    getline(*this->first_layer.get(), line); this->header.kmer_size = stoi(split(line, ' ')[1]);
    getline(*this->first_layer.get(), line); this->header.nb_counts = stoi(split(line, ' ')[1]);
    getline(*this->first_layer.get(), line); this->header.is_hashes = stoi(split(line, ' ')[1]);
    getline(*this->first_layer.get(), line); this->header.is_compressed = stoi(split(line, ' ')[1]);
    getline(*this->first_layer.get(), line); this->header.second_magic = stoi(split(line, ' ')[1]);
  }

public:
  CountMatrixFile() = delete;
  CountMatrixFile(CountMatrixFile const &) = delete;
  CountMatrixFile& operator= (CountMatrixFile const &) = delete;
  CountMatrixFile(CountMatrixFile&&) = delete;
  CountMatrixFile& operator= (CountMatrixFile&&) = delete;

  ~CountMatrixFile() = default;

  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, is>{}, void>::type>
  CountMatrixFile(const string& path)
    : IFile<cmheader_t, stream, buf_size>(
      path, cmheader_d, std::ios::binary | std::ios::in)
  {
    if (MatrixType != matrix_format_t::ASCII)
      this->first_layer->read(reinterpret_cast<char*>(&this->header), sizeof(this->header));
    else
      read_ascii_header();

    if (this->header.first_magic != magic1 || this->header.second_magic != magic2)
      throw std::runtime_error("Unable to read " + this->path + ". Possibly due to bad format.");
    
    std::stringstream msg;
    msg << "Unable to read "
        << fmt_to_string.at(this->header.matrix_type)
        << " matrix using MatrixFile<"
        << fmt_to_string.at(MatrixType) << ">.";
    if (MatrixType != this->header.matrix_type)
      throw std::runtime_error(msg.str());

    this->template set_second_layer<icstream>(this->header.is_compressed); 
  }

  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, os>{}, void>::type>
  CountMatrixFile(const string& path, int32_t partition_id,
                  uint32_t nb_counts, uint32_t kmer_size,
                  uint32_t is_hashes, uint32_t is_compressed)
    : IFile<cmheader_t, stream, buf_size>(path, cmheader_d, std::ios::binary | std::ios::out)
  {
    this->header.matrix_type = MatrixType;
    this->header.partition_id = partition_id;
    this->header.nb_counts = nb_counts;
    this->header.kmer_size = kmer_size;
    this->header.is_compressed = is_compressed;
    this->header.is_hashes = is_hashes;

    if (MatrixType == matrix_t::ASCII)
      write_ascii_header();
    else
      this->first_layer->write(reinterpret_cast<char*>(&this->header), sizeof(this->header));
    this->first_layer->flush();

    this->template set_second_layer<ocstream>(this->header.is_compressed);
  }

  template<typename S = stream,
           matrix_t M = MatrixType,
           typename = typename std::enable_if<std::is_same<S, os>{}, void>::type>
  void write(Kmer<KmerType>& kmer, std::vector<CountType>& counts)
  {
    if (MatrixType == matrix_t::BIN)
    {
      this->second_layer->write(
        reinterpret_cast<char*>(kmer.get_data()), sizeof(KmerType));
      this->second_layer->write(
        reinterpret_cast<char*>(counts.data()), counts.size()*sizeof(CountType));
    }
    else
    {
      if (this->header.is_hashes)
        *this->second_layer.get() << *kmer.get_data();
      else
        *this->second_layer.get() << kmer.str_value();
      for (auto& c: counts)
        *this->second_layer.get() << " " << std::to_string(c);
      *this->second_layer.get() << "\n";
    }
  }

  template<typename S = stream,
           matrix_t M = MatrixType,
           typename = typename std::enable_if<std::is_same<S, os>{}, void>::type>
  void write(KmerType& kmer, std::vector<CountType>& counts)
  {
    if (MatrixType == matrix_t::BIN)
    {
      this->second_layer->write(
        reinterpret_cast<char*>(&kmer), sizeof(KmerType));
      this->second_layer->write(
        reinterpret_cast<char*>(counts.data()), counts.size()*sizeof(CountType));
    }
    else
    {
      *this->second_layer.get() << kmer;
      for (auto& c: counts)
        *this->second_layer.get() << " " << to_string(c);
      *this->second_layer.get() << "\n";
    }
  }

  template<typename S = stream,
           matrix_t M = MatrixType,
           typename = typename std::enable_if<std::is_same<S, is>{}, void>::type>
  bool read(Kmer<KmerType>& kmer, std::vector<CountType>& counts)
  {
    if (MatrixType == matrix_t::BIN)
    {
      KmerType tmp;
      this->second_layer->read(
        reinterpret_cast<char*>(&tmp), sizeof(KmerType));
      if (this->second_layer->eof()) return false;
      this->second_layer->read(
        reinterpret_cast<char*>(counts.data()), counts.size()*sizeof(CountType));
      kmer.set_kmer(tmp, this->header.kmer_size);
      return true;
    }
    else
    {
      string line;
      if (getline(*this->second_layer.get(), line))
      {
        std::vector<std::string> tmp = std::move(split(line, ' '));
        kmer.set_kmer(tmp[0]);
        for (size_t i=1; i<counts.size()+1; i++)
          counts[i-1] = stoi(tmp[i]);
        return true;
      }
      return false;
    }
  }

  template<typename S = stream,
           matrix_t M = MatrixType,
           typename = typename std::enable_if<std::is_same<S, is>{}, void>::type>
  bool read(KmerType& kmer, std::vector<CountType>& counts)
  {
    if (MatrixType == matrix_t::BIN)
    {
      this->second_layer->read(
        reinterpret_cast<char*>(&kmer), sizeof(KmerType));
      if (this->second_layer->eof()) return false;
      this->second_layer->read(
        reinterpret_cast<char*>(counts.data()), counts.size()*sizeof(CountType));
      return true;
    }
    else
    {
      string line;
      if (getline(*this->second_layer.get(), line))
      {
        std::vector<std::string> tmp = std::move(split(line, ' '));
        kmer = stoi(tmp[0]);
        for (size_t i=1; i<counts.size()+1; i++)
          counts[i-1] = stoi(tmp[i]);
        return true;
      }
      return false;
    }
  }
};

//! \ingroup IO
//! \class PAMatrix
//! \brief TODO
template<typename stream,
         typename KmerType,
         size_t buf_size = 4096,
         typename = typename std::enable_if<std::is_integral<KmerType>::value ||
                                            std::is_same<KmerType, __uint128_t>{}, void>::type>
class PAMatrixFile : public IFile<paheader_t, stream, buf_size>
{

  using ocstream = lz4_stream::basic_ostream<buf_size>;
  using icstream = lz4_stream::basic_istream<buf_size>;

public:
  PAMatrixFile() = delete;
  PAMatrixFile(PAMatrixFile const &) = delete;
  PAMatrixFile& operator= (PAMatrixFile const &) = delete;
  PAMatrixFile(PAMatrixFile&&) = delete;
  PAMatrixFile& operator= (PAMatrixFile&&) = delete;

  ~PAMatrixFile() = default;

  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, is>{}, void>::type>
  PAMatrixFile(const string& path) 
    : IFile<paheader_t, stream, buf_size>(path, paheader_d, std::ios::binary | ios::in)
  {
    this->first_layer->read(reinterpret_cast<char*>(&this->header), sizeof(this->header));
    if (this->header.first_magic != magic1 || this->header.second_magic != magic2)
      throw std::runtime_error("Unable to read " + this->path + ". Possibly due to bad format.");

    if (this->header.ktsize != sizeof(KmerType))
      throw std::runtime_error("Invalid kmer type size");
    this->template set_second_layer<icstream>(this->header.is_compressed);
  }

  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, os>{}, void>::type>
  PAMatrixFile(const string& path, int32_t partition_id,
               uint32_t nb_files, uint32_t kmer_size,
               uint32_t is_hashes, uint32_t is_compressed)
    : IFile<paheader_t, stream, buf_size>(path, paheader_d, std::ios::binary | ios::out)
  {
    this->header.partition_id = partition_id;
    this->header.ktsize = sizeof(KmerType);
    this->header.bits_in_use = nb_files;
    this->header.size_in_bytes = NBYTE(nb_files);
    this->header.kmer_size = kmer_size;
    this->header.is_hashes = is_hashes;
    this->header.is_compressed = is_compressed;

    this->first_layer->write(reinterpret_cast<char*>(&this->header), sizeof(this->header)); 
    this->first_layer->flush();
    this->template set_second_layer<ocstream>(this->header.is_compressed);
  }

  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, os>{}, void>::type>
  void write(Kmer<KmerType>& kmer, std::vector<char>& bit_vector)
  {
    this->second_layer->write(reinterpret_cast<char*>(kmer.get_data()), sizeof(KmerType));
    this->second_layer->write(reinterpret_cast<char*>(bit_vector.data()), bit_vector.size());
  }
  
  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, os>{}, void>::type>
  void write(KmerType& kmer, std::vector<char>& bit_vector)
  {
    this->second_layer->write(reinterpret_cast<char*>(&kmer), sizeof(KmerType));
    this->second_layer->write(reinterpret_cast<char*>(bit_vector.data()), bit_vector.size());
  }

  template<typename data,
           typename S = stream,
           typename = typename std::enable_if<std::is_same<S, os>{}, void>::type>
  void write(KmerType& kmer, data* bit_vector, size_t len)
  {
    this->second_layer->write(reinterpret_cast<char*>(&kmer), sizeof(KmerType));
    this->second_layer->write(reinterpret_cast<char*>(bit_vector), len*sizeof(data));
  }

  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, is>{}, void>::type>
  bool read(Kmer<KmerType>& kmer, std::vector<char>& bit_vector)
  {
    this->second_layer->read(reinterpret_cast<char*>(&tmp), sizeof(KmerType));
    if (this->second_layer->eof()) return false;
    this->second_layer->read(reinterpret_cast<char*>(bit_vector.data()), bit_vector.size());
    kmer.set_kmer(tmp, this->header.kmer_size);
    return true;
  }

  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, is>{}, void>::type>
  bool read(KmerType& kmer, std::vector<char>& bit_vector)
  {
    this->second_layer->read(reinterpret_cast<char*>(&kmer), sizeof(KmerType));
    if (this->second_layer->eof()) return false;
    this->second_layer->read(reinterpret_cast<char*>(bit_vector.data()), bit_vector.size());
    return true;
  }

private:
  KmerType   tmp    {0};
};

template<typename stream,
         matrix_t MatrixType,
         size_t   buf_size = 4096,
         typename = typename std::enable_if<MatrixType == matrix_t::BF ||
                                            MatrixType == matrix_t::BIT, void>::type>
class BitMatrixFile : public IFile<bmheader_t, stream, buf_size>
{
  using ocstream  = lz4_stream::basic_ostream<buf_size>;
  using icstream  = lz4_stream::basic_istream<buf_size>;
public:
  BitMatrixFile() = delete;
  BitMatrixFile(BitMatrixFile const &) = delete;
  BitMatrixFile& operator= (BitMatrixFile const &) = delete;
  BitMatrixFile(BitMatrixFile&&) = delete;
  BitMatrixFile& operator= (BitMatrixFile&&) = delete;

  ~BitMatrixFile() = default;

  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, is>{}, void>::type>
  BitMatrixFile(const string& path)
    : IFile<bmheader_t, stream, buf_size>(
        path, bmheader_d, std::ios::binary | std::ios::in)
  {
    this->first_layer->read(reinterpret_cast<char*>(&this->header), sizeof(this->header));
    if (this->header.first_magic != magic1 || this->header.second_magic != magic2)
      throw std::runtime_error("Unable to read " + this->path + ". Possibly due to bad format.");

    //if (this->header.matrix_type != MatrixType)
    //  if (this->header.matrix_type != matrix_t::BF && MatrixType != matrix_t::BIT)
    //    throw std::runtime_error("Bad matrix type.");

    row_count = this->header.nb_rows_use;
    this->template set_second_layer<icstream>(this->header.is_compressed);
  }

  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, os>{}, void>::type>
  BitMatrixFile(const string& path, uint32_t partition_id,
                uint64_t nb_rows, uint64_t nb_cols,
                uint64_t is_compressed)
    : IFile<bmheader_t, stream, buf_size>(
        path, bmheader_d, std::ios::binary | std::ios::out)
  {
    this->header.matrix_type = MatrixType;
    this->header.partition_id = partition_id;
    this->header.nb_rows = nb_rows;
    this->header.nb_rows_use = nb_rows % 8 ? NMOD8(nb_rows) : nb_rows;
    this->header.nb_cols = nb_cols;
    this->header.nb_cols_use = nb_cols % 8 ? NMOD8(nb_cols) : nb_cols;
    this->header.size_in_bytes = NBYTE(this->header.nb_cols_use)*this->header.nb_rows_use;

    this->first_layer->write(reinterpret_cast<char*>(&this->header), sizeof(this->header)); 
    this->first_layer->flush();
    this->template set_second_layer<ocstream>(this->header.is_compressed);
  }

  template<matrix_t fmt = MatrixType,
           typename S = stream,
           typename = typename std::enable_if<fmt == matrix_format_t::BIT, void>::type,
           typename = typename std::enable_if<std::is_same<S, os>{}, S>::type>
  void dump(BitMatrix& bitmatrix)
  {
    if (this->header.size_in_bytes != bitmatrix.get_size_in_byte())
      throw std::runtime_error("Invalid BitMatrix size");
    this->second_layer->write(
      reinterpret_cast<char*>(bitmatrix.matrix), bitmatrix.get_size_in_byte()
    );
  }

  template<matrix_t fmt = MatrixType,
           typename S = stream,
           typename = typename std::enable_if<fmt == matrix_format_t::BIT, void>::type,
           typename = typename std::enable_if<std::is_same<S, is>{}, S>::type>
  void load(BitMatrix& bitmatrix)
  {
    if (this->header.size_in_bytes != bitmatrix.get_size_in_byte())
      throw std::runtime_error("Invalid BitMatrix size");
    this->second_layer->read(
      reinterpret_cast<char*>(bitmatrix.matrix), bitmatrix.get_size_in_byte()
    );
  }

  template<matrix_t fmt = MatrixType,
           typename S = stream,
           typename = typename std::enable_if<fmt == matrix_format_t::BF, void>::type,
           typename = typename std::enable_if<std::is_same<S, os>{}, S>::type>
  bool write(vector<char>& bit_vector)
  {
    if (row_count >= this->header.nb_rows)
      return false;
    this->second_layer->write(bit_vector.data(), bit_vector.size());
    row_count++;
    return true;
  }

  template<typename data,
           matrix_t fmt = MatrixType,
           typename S = stream,
           typename = typename std::enable_if<fmt == matrix_format_t::BF, void>::type,
           typename = typename std::enable_if<std::is_same<S, os>{}, S>::type>
  bool write(data* bit_vector, size_t len)
  {
    if (row_count >= this->header.nb_rows)
      return false;
    this->second_layer->write(reinterpret_cast<char*>(bit_vector), len*sizeof(data));
    row_count++;
    return true;
  }
  
  template<matrix_t fmt = MatrixType,
           typename S = stream,
           typename = typename std::enable_if<fmt == matrix_format_t::BF, void>::type,
           typename = typename std::enable_if<std::is_same<S, is>{}, S>::type>
  bool read(std::vector<char>& bit_vector)
  {
    this->second_layer->read(
      reinterpret_cast<char*>(bit_vector.data()), bit_vector.size());
    if (this->second_layer->eof()) return false;
    return true;
  }

  template<typename data,
           matrix_t fmt = MatrixType,
           typename S = stream,
           typename = typename std::enable_if<fmt == matrix_format_t::BF, void>::type,
           typename = typename std::enable_if<std::is_same<S, is>{}, S>::type>
  bool read(data* bit_vector, size_t len)
  {
    this->second_layer->read(
      reinterpret_cast<char*>(bit_vector), sizeof(data)*len);
    if (this->second_layer->eof()) return false;
    return true;
  }

  template<matrix_t fmt = MatrixType,
           typename = typename std::enable_if<fmt == matrix_t::BF, void>::type>
  bool is_consistent()
  {
    return row_count == this->header.nb_rows_use;
  }

  uint64_t get_row_count()
  {
    return row_count;
  }

private:
  uint64_t row_count {0};
};

//! \ingroup IO
//! \class BitVectorFile
//! \brief TODO
template<typename stream   = os,
         size_t   buf_size = 4096>
class BitVectorFile : public IFile<bvheader_t, stream, buf_size>
{
  using ocstream = lz4_stream::basic_ostream<buf_size>;
  using icstream = lz4_stream::basic_istream<buf_size>;
  using buffer_t = std::vector<char>;

public:
  BitVectorFile () = delete;
  BitVectorFile (BitVectorFile const &) = delete;
  BitVectorFile& operator= (BitVectorFile const &) = delete;
  BitVectorFile (BitVectorFile&&) = delete;
  BitVectorFile& operator= (BitVectorFile&&) = delete;

  ~BitVectorFile() = default;

  //! \brief
  //! \param :
  //! \param :
  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, is>{}, S>::type>
  BitVectorFile(const std::string& path)
    : IFile<bvheader_t, stream, buf_size>(
        path, bvheader_d, std::ios::binary | std::ios::in)
  {
    this->first_layer->read(reinterpret_cast<char*>(&this->header), sizeof(this->header));
    if (this->header.first_magic != magic1
        || this->header.second_magic != magic2)
      throw std::runtime_error("Unable to read: " + path + ". Possibly due to bad format.");

    this->template set_second_layer<icstream>(this->header.is_compressed); 
  }

  //! \brief
  //! \param :
  //! \param :
  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, os>{}, S>::type>
  BitVectorFile(const std::string& path, uint32_t file_id,
                uint32_t partition_id, uint64_t nb_bits, 
                uint32_t is_compressed) 
    : IFile<bvheader_t, stream, buf_size>(path, bvheader_d, std::ios::binary | std::ios::out)
  {
    if (nb_bits % 8 != 0)
      throw runtime_error("nb_bits must be a multiple of 8");
    this->header.file_id        = file_id;
    this->header.partition_id   = partition_id;
    this->header.partition_size = NBYTE(nb_bits);
    this->header.nb_bits        = nb_bits;
    this->header.is_compressed  = is_compressed;

    this->first_layer->write(reinterpret_cast<char*>(&this->header), sizeof(this->header));
    this->first_layer->flush();

    this->template set_second_layer<ocstream>(this->header.is_compressed);
  }

  //! \brief
  //! \param :
  //! \param :
  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, os>{}, S>::type>
  bool write(std::vector<char>& bit_vector)
  {
    if (!is_rw)
    {
      this->second_layer->write(reinterpret_cast<char*>(bit_vector.data()), bit_vector.size());
      is_rw = true;
      return true;
    }
    return false;
  }

  //! \brief
  //! \param :
  //! \param :
  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, is>{}, S>::type>
  bool read(std::vector<char>& bit_vector)
  {
    if (bit_vector.size() != this->header.partition_size)
      throw std::runtime_error("Provided bit-vector size and partition size differ.");
    if (!is_rw)
    {
      this->second_layer->read(reinterpret_cast<char*>(bit_vector.data()), bit_vector.size());
      is_rw = true;
    }
    return false;
  }

  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, is>{}, S>::type>
  bool read(char* bit_vector, size_t len)
  {
    if (len != this->header.partition_size)
      throw std::runtime_error("Provided bit-vector size and partition size differ.");
    if (!is_rw)
    {
      this->second_layer->read(bit_vector, len);
      is_rw = true;
    }
    return false;
  }

  //! \brief
  //! \param :
  //! \param :
  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, is>{}, void>::type>
  std::vector<char> read()
  {
    if (!is_rw)
    {
      std::vector<char> v(this->header.partition_size);
      this->second_layer->read(reinterpret_cast<char*>(v.data()), v.size());
      is_rw = true;
      return v;
    }
    return vector<char>(0);
  }

  //! \brief
  //! \param :
  //! \param :
  std::pair<uint64_t, uint64_t> get_window()
  {
    uint64_t lower = this->header.partition_size*this->header.partition_id*8;
    uint64_t upper = lower + this->header.partition_size*8;
    return {lower, upper-1};
  }

private:
  bool is_rw  {false};
};

template<typename stream,
         size_t buf_size = 4096>
class HistFile : public IFile<hheader_t, stream, buf_size>
{
  using ocstream = lz4_stream::basic_ostream<buf_size>;
  using icstream = lz4_stream::basic_istream<buf_size>;
  using buffer_t = std::vector<char>;

public:
  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, is>{}, S>::type>
  HistFile(const std::string& path)
    : IFile<hheader_t, stream, buf_size>(
        path, hheader_d, std::ios::binary | std::ios::in)
  {
    this->first_layer->read(reinterpret_cast<char*>(&this->header), sizeof(this->header));
    if (this->header.first_magic != magic1
        || this->header.second_magic != magic2)
      throw std::runtime_error("Unable to read: " + path + ". Possibly due to bad format.");

    this->template set_second_layer<icstream>(false); 
  }

  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, os>{}, S>::type>
  HistFile(const KHist& hist, const std::string& path)
    : IFile<hheader_t, stream, buf_size>(path, hheader_d, std::ios::binary | std::ios::out)
  {
    this->header.id        = hist.idx;
    this->header.kmer_size = hist.ksize;
    this->header.lower     = hist.lower;
    this->header.upper     = hist.upper;
    this->header.uniq      = hist.uniq;
    this->header.total     = hist.total;

    this->first_layer->write(reinterpret_cast<const char*>(&this->header), sizeof(this->header));
    this->first_layer->flush();

    this->template set_second_layer<ocstream>(false);

    this->second_layer->write(reinterpret_cast<const char*>(&hist.oob_lu), sizeof(uint64_t));
    this->second_layer->write(reinterpret_cast<const char*>(&hist.oob_uu), sizeof(uint64_t));
    this->second_layer->write(reinterpret_cast<const char*>(&hist.oob_ln), sizeof(uint64_t));
    this->second_layer->write(reinterpret_cast<const char*>(&hist.oob_un), sizeof(uint64_t));
    
    this->second_layer->write(
      reinterpret_cast<const char*>(hist.hist_u.data()), hist.hist_u.size()*sizeof(uint64_t));
    this->second_layer->write(
      reinterpret_cast<const char*>(hist.hist_n.data()), hist.hist_n.size()*sizeof(uint64_t));
  }

  template<typename S = stream,
           typename = typename std::enable_if<std::is_same<S, is>{}, void>::type>
  KHist read()
  {
    KHist hist(this->header.id, this->header.kmer_size, this->header.lower, this->header.upper);
    
    this->second_layer->read(reinterpret_cast<char*>(&hist.oob_lu), sizeof(uint64_t));
    this->second_layer->read(reinterpret_cast<char*>(&hist.oob_uu), sizeof(uint64_t));
    this->second_layer->read(reinterpret_cast<char*>(&hist.oob_ln), sizeof(uint64_t));
    this->second_layer->read(reinterpret_cast<char*>(&hist.oob_un), sizeof(uint64_t));

    this->second_layer->read(
      reinterpret_cast<char*>(hist.hist_u.data()), hist.hist_u.size()*sizeof(uint64_t));
    this->second_layer->read(
      reinterpret_cast<char*>(hist.hist_n.data()), hist.hist_n.size()*sizeof(uint64_t));
    return hist;
  }
};
} // end of namespace km
