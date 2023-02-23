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
#include <vector>
#include <filesystem>
#include <memory>
#include <kmtricks/io/superk_file.hpp>
#include <gatb/gatb_core.hpp>
#include <gatb/system/api/IThread.hpp>
#include <unordered_set>

namespace fs = std::filesystem;

namespace km {

class SuperKStorageReader
{
public:
  SuperKStorageReader() {}
  SuperKStorageReader(const std::string& prefix)
  {
    std::ifstream info(fmt::format("{}/SuperKmerBinInfoFile", prefix));
    check_fstream_good(fmt::format("{}/SuperKmerBinInfoFile", prefix), info);
    std::string line;

    std::getline(info, line); m_base = line;
    std::getline(info, line); m_path = line;
    std::getline(info, line); m_nb_files = std::stoi(line);

    m_nbk_per_file.resize(m_nb_files, 0);
    m_file_size.resize(m_nb_files, 0);
    m_files.resize(m_nb_files, nullptr);
    m_synchros.resize(m_nb_files, 0);

    for (size_t i=0; i<m_files.size(); i++)
    {
      std::getline(info, line); m_nbk_per_file[i] = std::stoll(line);
      std::getline(info, line); m_file_size[i] = std::stoll(line);
    }
  }

  ~SuperKStorageReader()
  {
    closeFiles();
  }

  void flushFile(int fileId)
  {
    m_files[fileId]->flush();
  }

  void flushFiles()
  {
    for (int i=0; i<m_nb_files; i++)
      flushFile(i);
  }

  void eraseFiles()
  {
    for (int i=0; i<m_nb_files; i++)
      eraseFile(i);
  }

  void eraseFile(int fileId)
  {
    closeFile(fileId);
    std::remove(fmt::format("{}/{}.{}", m_path, m_base, fileId).c_str());
  }

  void openFiles()
  {
    fs::create_directory(m_path);
    for (int i=0; i<m_nb_files; i++)
      openFile(i);
  }

  void openFile(int fileId)
  {
    std::string path = fmt::format("{}/{}.{}", m_path, m_base, fileId);
    m_files[fileId] = std::make_shared<SuperkReader<8192>>(path);
    m_synchros[fileId] = gatb::core::system::impl::System::thread().newSynchronizer();
    m_synchros[fileId]->use();
  }

  void closeFiles()
  {
    for (int i=0; i<m_nb_files; i++)
      closeFile(i);
  }

  void closeFile(int fileId)
  {
    if (!m_files.empty())
    {
      if (m_files[fileId])
      {
        m_files[fileId]->close();
        m_files[fileId] = nullptr;
        m_synchros[fileId]->forget();
      }
    }
  }

  int readBlock(unsigned char** block,
                unsigned int* max_block_size,
                unsigned int* nb_bytes_read,
                int file_id)
  {
    m_synchros[file_id]->lock();
    int nbr = m_files[file_id]->read_size(nb_bytes_read);

    if (nbr == 0)
    {
      m_synchros[file_id]->unlock();
      return 0;
    }
    if (*nb_bytes_read > *max_block_size)
    {
      *block = (unsigned char*) realloc(*block, *nb_bytes_read);
      *max_block_size = *nb_bytes_read;
    }
    m_files[file_id]->read_block(*block, *nb_bytes_read);

    m_synchros[file_id]->unlock();
    return *nb_bytes_read;
  }

  int nbFiles() const { return m_nb_files; }

  std::string getFileName(int fileId) const
  {
    return fmt::format("{}/{}.{}", m_path, m_base, fileId);
  }

  uint64_t getNbItems(int fileId) const
  {
    return m_nbk_per_file[fileId];
  }

  uint64_t getFileSize(int fileId) const
  {
    return m_file_size[fileId];
  }

private:
  std::string m_base;
  std::string m_path;
  std::vector<uint64_t> m_nbk_per_file;
  std::vector<uint64_t> m_file_size;
  std::vector<skr_t<8192>> m_files;
  std::vector<gatb::core::system::ISynchronizer*> m_synchros;
  int m_nb_files;
};

using sk_storage_t = std::shared_ptr<SuperKStorageReader>;

class SuperKStorageWriter
{
public:
  SuperKStorageWriter(const std::string& prefix, const std::string& name,
                      size_t nb_files, bool lz4, std::unordered_set<int> restricted)
    : m_path(prefix), m_base(name), m_nb_files(nb_files), m_lz4(lz4), m_restricted(restricted)
  {
    m_nbk_per_file.resize(m_nb_files, 0);
    m_file_size.resize(m_nb_files, 0);
    m_files.resize(m_nb_files, nullptr);
    m_synchros.resize(m_nb_files, nullptr);
    openFiles();

    m_capacity = 32768;
    m_max_superk_size = 255;
    m_buffers.resize(m_nb_files);
    m_buffers_idx.resize(m_nb_files, 0);
    for (unsigned int ii=0; ii<m_buffers.size(); ii++)
    {
      m_buffers[ii] = reinterpret_cast<uint8_t*>(MALLOC(sizeof(uint8_t) * m_capacity));
    }
  }

  void flushAllCache()
  {
    for (unsigned int ii=0; ii<m_buffers.size(); ii++)
    {
      flushCache(ii);
    }
  }

  void flushCache(int file_id)
  {
    if (m_buffers_idx[file_id] != 0)
    {
      writeBlock(m_buffers[file_id], m_buffers_idx[file_id], file_id, m_nbk_per_file[file_id]);
      m_buffers_idx[file_id] = 0;
      m_nbk_per_file[file_id] = 0;
    }
  }

  void insertSuperkmer(uint8_t* superk, int nb_bytes, uint8_t nbk, int file_id)
  {
    if ( (m_buffers_idx[file_id]+nb_bytes+1) > static_cast<int>(m_capacity))
    {
      flushCache(file_id);
    }
    m_buffers[file_id][m_buffers_idx[file_id]++] = nbk;
    memcpy(m_buffers[file_id] + m_buffers_idx[file_id], superk, nb_bytes);
    m_buffers_idx[file_id] += nb_bytes;
    m_nbk_per_file[file_id] += nbk;
  }

  ~SuperKStorageWriter()
  {
    flushAllCache();
    for (unsigned int ii=0; ii<m_buffers.size(); ii++)
      FREE(m_buffers[ii]);
    closeFiles();
  }

  void flushFile(int fileId)
  {
    m_files[fileId]->flush();
  }

  void flushFiles()
  {
    for (int i=0; i<m_nb_files; i++)
      flushFile(i);
  }

  void eraseFiles()
  {
    for (int i=0; i<m_nb_files; i++)
      eraseFile(i);
  }

  void eraseFile(int fileId)
  {
    closeFile(fileId);
    std::remove(fmt::format("{}/{}.{}", m_path, m_base, fileId).c_str());
  }

  void openFiles()
  {
    fs::create_directory(m_path);
    for (int i=0; i<m_nb_files; i++)
      openFile(i);
  }

  void openFile(int fileId)
  {
    if (!m_restricted.count(fileId))
      return;

    std::string path = fmt::format("{}/{}.{}", m_path, m_base, fileId);
    m_files[fileId] = std::make_shared<SuperkWriter<8192>>(path,
                                                           fileId,
                                                           m_lz4);
    m_synchros[fileId] = gatb::core::system::impl::System::thread().newSynchronizer();
    m_synchros[fileId]->use();
  }

  void closeFiles()
  {
    for (int i=0; i<m_nb_files; i++)
      closeFile(i);
  }

  void closeFile(int fileId)
  {
    if (m_files[fileId])
    {
      if (!m_restricted.count(fileId))
        return;
      m_files[fileId]->close();
      m_files[fileId] = nullptr;
      m_synchros[fileId]->forget();
    }
  }

  void writeBlock(unsigned char* block,
                unsigned int block_size,
                int file_id,
                int nbkmers)
  {
    if (!m_restricted.count(file_id))
      return;
    m_synchros[file_id]->lock();
    m_nbk_per_file[file_id] += nbkmers;
    m_file_size[file_id] += block_size + sizeof(block_size);
    m_files[file_id]->write_size(block_size);
    m_files[file_id]->write_block(block, block_size);
    m_synchros[file_id]->unlock();
  }

  int nbFiles() const { return m_nb_files; }

  std::string getFileName(int fileId) const
  {
    return fmt::format("{}.{}", m_base, fileId);
  }

  uint64_t getNbItems(int fileId) const
  {
    return m_nbk_per_file[fileId];
  }

  uint64_t getFileSize(int fileId) const
  {
    return m_file_size[fileId];
  }

  void SaveInfoFile(const std::string& prefix)
  {
    std::ofstream info(fmt::format("{}/SuperKmerBinInfoFile", prefix));
    check_fstream_good(fmt::format("{}/SuperKmerBinInfoFile", prefix), info);
    info << m_base << "\n";
    info << m_path << "\n";
    info << m_nb_files << "\n";
    for (int i=0; i<m_nb_files; i++)
    {
      info << m_nbk_per_file[i] << "\n";
      info << m_file_size[i] << "\n";
    }
  }
private:
  std::string m_base;
  std::string m_path;
  std::vector<uint64_t> m_nbk_per_file;
  std::vector<uint64_t> m_file_size;
  std::vector<skw_t<8192>> m_files;
  std::vector<gatb::core::system::ISynchronizer*> m_synchros;
  std::unordered_set<int> m_restricted;
  int m_nb_files;
  bool m_lz4;

  size_t m_max_superk_size;
  size_t m_capacity;
  std::vector<uint8_t*> m_buffers;
  std::vector<int> m_buffers_idx;
};

};
