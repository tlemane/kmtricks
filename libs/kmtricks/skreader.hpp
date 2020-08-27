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
#include <fstream>
#include <iostream>
#include "sequences.hpp"

using namespace std;

namespace km
{

class SuperkStorage
{
public:
  SuperkStorage(string &superk_dir, string &part_prefix, int nb_parts);
  ~SuperkStorage();

  void reset_all();
  void reset(int part_id);
  void close_files();
  void close_file(int part_id);

  int read_block(uchar **block, uint *max_block_size, uint *nb_bytes, int id);

  int nb_files();

public:
  vector<ifstream> _parts;
private:
  string _pdir;
  int    _nb_parts;
  string _prefix;
};

#ifndef _KM_LIB_INCLUDE_
SuperkStorage::SuperkStorage(string &superk_dir, string &part_prefix, int nb_parts)
: _pdir(superk_dir),
  _nb_parts(nb_parts),
  _prefix(part_prefix)
{
  string path = _pdir + "/" + _prefix;
  for (int i=0; i<_nb_parts; i++)
  {
    string ppath = path + to_string(i);
    _parts.push_back(ifstream(ppath, ios::in | ios::binary));
  }
}


SuperkStorage::~SuperkStorage()
{
  close_files();
}


void SuperkStorage::reset_all()
{
  for (auto & file : _parts)
  {
    file.clear();
    file.seekg(0, ios::beg);
  }
}


void SuperkStorage::reset(int part_id)
{
  _parts[part_id].clear();
  _parts[part_id].seekg(0, ios::beg);
}


void SuperkStorage::close_files()
{
  for (auto & file : _parts)
    file.close();
}


void SuperkStorage::close_file(int part_id)
{
  _parts[part_id].close();
}


int SuperkStorage::read_block(uchar **block, uint *max_block_size, uint *nb_bytes, int id)
{
  _parts[id].read((char*)nb_bytes, sizeof(*max_block_size));
  int nb = _parts[id].gcount();

  if (nb == 0)
  {
    *nb_bytes = 0;
    return 0;
  }

  if (*nb_bytes > *max_block_size)
  {
    *block = static_cast<uchar*>(realloc(*block, *nb_bytes));
    *max_block_size = *nb_bytes;
  }

  _parts[id].read((char*)*block, *nb_bytes);

  return *nb_bytes;
}


int SuperkStorage::nb_files()
{
  return _nb_parts;
}
#endif

template<typename K>
class SuperkReader
{
public:
  SuperkReader(SuperkStorage *sk_storage, size_t kmer_size);
  ~SuperkReader();

  bool next_superk(int part_id, Superk<K> *superk);

private:
  size_t          _ksize;
  SuperkStorage   *_sk_storage;
  vector<uchar*>  _buffers;
  vector<uint>    _buffer_sizes;
  vector<uint>    _nb_bytes_reads;
  vector<uchar*>  _current;
  vector<bool>    _init;
};

#ifndef _KM_LIB_INCLUDE_

template<typename K>
SuperkReader<K>::SuperkReader(SuperkStorage *sk_storage, size_t kmer_size)
: _ksize(kmer_size),
  _sk_storage(sk_storage)
{
  int nb_files = _sk_storage->nb_files();
  _buffers.resize(nb_files);
  _buffer_sizes.resize(nb_files);
  _nb_bytes_reads.resize(nb_files);
  _current.resize(nb_files);
  _init.resize(nb_files, false);
}


template<typename K>
SuperkReader<K>::~SuperkReader<K>()
{
  for (auto &buf : _buffers)
    if (buf)
      free(buf);
}


template<typename K>
bool SuperkReader<K>::next_superk(int part_id, Superk<K> *superk)
{
  if (!_init[part_id])
  {
    _sk_storage->read_block(&_buffers[part_id], &_buffer_sizes[part_id], &_nb_bytes_reads[part_id], part_id);
    _init[part_id] = true;
    _current[part_id] = _buffers[part_id];
  }
  if (!(_current[part_id] < (_buffers[part_id]+_nb_bytes_reads[part_id])))
  {
    _sk_storage->read_block(&_buffers[part_id], &_buffer_sizes[part_id], &_nb_bytes_reads[part_id], part_id);
    _current[part_id] = _buffers[part_id];
  }

  if (!_nb_bytes_reads[part_id])
  {
    free(_buffers[part_id]);
    _buffers[part_id] = nullptr;
    return false;
  }

  uint8_t nbk = *_current[part_id]; _current[part_id]++;
  size_t superk_size = nbk + _ksize - 1;
  superk->set_superk(_current[part_id], superk_size, _ksize, true);

  _current[part_id] += (_ksize+nbk-1+3)/4;
  return true;
}
#endif
} // end of namespace km