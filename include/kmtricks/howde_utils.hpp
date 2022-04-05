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
#include <sdsl/bit_vectors.hpp>
#include <bloom_filter_file.h>
#include <kmtricks/kmdir.hpp>
#include <kmtricks/io/vector_file.hpp>
#include <kmtricks/io/vector_matrix_file.hpp>
#include <kmtricks/hash.hpp>

#define _FILE_OFFSET_BITS 64

#ifdef __linux__
  #include <linux/version.h>
  #if LINUX_VERSION_CODE >= KERNEL_VERSION(4, 5, 0)
    #include <cfrcat/cfrcat.hpp>
  #else
    #include <sys/sendfile.h>
  #endif
#elif __APPLE__
  #include <sys/types.h>
  #include <sys/socket.h>
  #include <sys/uio.h>
#endif

#define round_up_16(b)  ((((std::uint64_t) (b))+15)&(~15))

namespace km {

class IBloomBuilder
{
public:
  IBloomBuilder(OUT_FORMAT bf_type, uint64_t bloom_size, uint32_t file_id, uint32_t nb_parts, uint32_t kmer_size)
    : m_bf_type(bf_type), m_bloom_size(bloom_size), m_file_id(file_id), m_kmer_size(kmer_size), m_nb_parts(nb_parts)
  {
    m_hw = HashWindow(KmDir::get().m_hash_win);
  }

protected:
  void write_header(std::ostream& stream)
  {
    uint32_t header_size = round_up_16(bffileheader_size(1));
    bffileheader* header = reinterpret_cast<bffileheader*>(new char[header_size]());
    header->magic = bffileheaderMagicUn;
    header->headerSize = sizeof(bffileprefix);
    stream.write(reinterpret_cast<char*>(header), header_size);
    size_t bw = header_size;

    header->magic = bffileheaderMagic;
    header->headerSize = header_size;
    header->version = bffileheaderVersion;
    header->bfKind = bfkind_simple;
    header->smerSize = m_kmer_size;
    header->numHashes = 1;
    header->hashSeed1 = 0;
    header->hashSeed2 = 0;
    header->hashModulus = m_bloom_size;
    header->numBits = m_bloom_size;
    header->numVectors = 1;
    header->setSizeKnown = false;
    header->setSize = 0;

    header->info[0].compressor = bvcomp_uncompressed;
    header->info[0].name = 0;
    header->info[0].offset = bw;

    header->info[0].numBytes = (m_bloom_size / 8) + sizeof(uint64_t);
    header->info[0].filterInfo = (uint64_t)0;
    stream.seekp(std::ios::beg);
    stream.write(reinterpret_cast<char*>(header), header_size);
    delete[] header;
  }

  void write_header_fd(int fd)
  {
    uint32_t header_size = round_up_16(bffileheader_size(1));
    bffileheader* header = reinterpret_cast<bffileheader*>(new char[header_size]());
    header->magic = bffileheaderMagicUn;
    header->headerSize = sizeof(bffileprefix);
    write(fd, reinterpret_cast<char*>(header), header_size);
    size_t bw = header_size;

    header->magic = bffileheaderMagic;
    header->headerSize = header_size;
    header->version = bffileheaderVersion;
    header->bfKind = bfkind_simple;
    header->smerSize = m_kmer_size;
    header->numHashes = 1;
    header->hashSeed1 = 0;
    header->hashSeed2 = 0;
    header->hashModulus = m_bloom_size;
    header->numBits = m_bloom_size;
    header->numVectors = 1;
    header->setSizeKnown = false;
    header->setSize = 0;

    header->info[0].compressor = bvcomp_uncompressed;
    header->info[0].name = 0;
    header->info[0].offset = bw;

    header->info[0].numBytes = (m_bloom_size / 8) + sizeof(uint64_t);
    header->info[0].filterInfo = (uint64_t)0;
    lseek(fd, 0, SEEK_SET);
    write(fd, reinterpret_cast<char*>(header), header_size);
    delete[] header;
  }

protected:
  OUT_FORMAT m_bf_type;
  uint64_t m_bloom_size;
  HashWindow m_hw;
  uint32_t m_file_id;
  uint32_t m_nb_parts;
  uint32_t m_kmer_size;
};

class BloomBuilderFromHash : public IBloomBuilder
{
public:
  BloomBuilderFromHash(
    std::vector<int>& files, std::vector<std::mutex>& file_mutex,
    OUT_FORMAT bf_type, uint64_t bloom_size, uint32_t file_id, uint32_t nb_parts, uint32_t kmer_size)
    : IBloomBuilder(bf_type, bloom_size, file_id, nb_parts, kmer_size), m_fds(files), m_mutex(file_mutex)
  {
  }

  void build()
  {
    std::string out_path = KmDir::get().get_filter_path(KmDir::get().m_fof.get_id(this->m_file_id), this->m_bf_type);
    int out_fd = open(out_path.c_str(), O_CREAT | O_TRUNC | O_RDWR, 0x01B6);
    this->write_header_fd(out_fd);
    write(out_fd, reinterpret_cast<char*>(&this->m_bloom_size), sizeof(this->m_bloom_size));
    for (size_t p=0; p<this->m_nb_parts; p++)
    {
      {
        std::unique_lock<std::mutex> lock(m_mutex[p]);

        lseek(m_fds[p], 49 + m_file_id * m_hw.get_window_size_bytes(), SEEK_SET);

#ifdef __linux__
  #if LINUX_VERSION_CODE >= KERNEL_VERSION(4, 5, 0)
        copy_file_range(m_fds[p], nullptr, out_fd, nullptr, m_hw.get_window_size_bytes(), 0);
  #else
        sendfile(out_fd, m_fds[p], nullptr, m_hw.get_window_size_bytes());
  #endif
#elif __APPLE__
        off_t len = m_hw.get_window_size_bytes();
        uint64_t n = len / 8192;
        uint64_t r = len % 8192;

        char buffer[8192];
        for (uint64_t i=0; i<n; ++i)
        {
          read(m_fds[p], buffer, 8192);
          write(out_fd, buffer, 8192);
        }

        read(m_fds[p], buffer, r);
        write(out_fd, buffer, r);
#endif
      }
    }
    close(out_fd);
  }

private:
  std::vector<int> m_fds;
  std::vector<std::mutex>& m_mutex;
};

class BloomBuilderFromVec : public IBloomBuilder
{
public:
  BloomBuilderFromVec(uint32_t file_id, OUT_FORMAT bf_type, uint64_t bloom_size,
                       uint32_t nb_parts, uint32_t kmer_size, bool lz4)
    : IBloomBuilder(bf_type, bloom_size, file_id, nb_parts, kmer_size), m_lz4(lz4)
  {
  }

  void build()
  {
    std::string out_path = KmDir::get().get_filter_path(KmDir::get().m_fof.get_id(this->m_file_id), this->m_bf_type);
    std::ofstream out(out_path, std::ios::binary|std::ios::out); check_fstream_good(out_path, out);

    this->write_header(out);
    out.write(reinterpret_cast<char*>(&this->m_bloom_size), sizeof(this->m_bloom_size));

    std::vector<char> buffer(this->m_hw.get_window_size_bytes());
    for (size_t p=0; p<this->m_nb_parts; p++)
    {
      BitVectorReader bvr(KmDir::get().get_count_part_path(KmDir::get().m_fof.get_id(this->m_file_id), p, m_lz4, KM_FILE::VECTOR));
      bvr.read(&buffer[0], this->m_hw.get_window_size_bytes());
      out.write(&buffer[0], this->m_hw.get_window_size_bytes());
    }
  }
private:
  bool m_lz4;
};

};
