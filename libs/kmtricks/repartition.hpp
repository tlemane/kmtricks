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
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>

using namespace std;

namespace km
{


class RepartFile
{
public:
  RepartFile(string m_path, string f_path = "");
  ~RepartFile();

  RepartFile(const RepartFile &r);
  RepartFile& operator=(const RepartFile &r);

  void load();
  uint16_t get(uint64_t minimizer_value);

private:
  string            _path;
  string            _path_freq;
  uint16_t          _nb_part;
  uint64_t          _nb_minims;
  uint16_t          _nb_pass;
  bool              _has_minim_freq;
  uint32_t          _magic;
  uint32_t          *_freq_order;
  vector<uint16_t>  _repart_table;

public:
  bool is_load;
};

#ifdef _KM_LIB_INCLUDE_
extern const uint32_t MAGIC_NUMBER;
#endif

#ifndef _KM_LIB_INCLUDE_
const uint32_t MAGIC_NUMBER = 0x12345678;

RepartFile::RepartFile(string m_path, string f_path)
: _path(m_path),
  _path_freq(f_path),
  _nb_part(0),
  _nb_minims(0),
  _nb_pass(0),
  _has_minim_freq(false),
  _magic(0),
  _freq_order(nullptr),
  is_load(false)
{
  load();
}


RepartFile::~RepartFile()
{
  if (_freq_order) delete[] _freq_order;
}


RepartFile::RepartFile(const RepartFile &r)
: _path(r._path),
  _path_freq(r._path_freq),
  _nb_part(r._nb_part),
  _nb_minims(r._nb_minims),
  _has_minim_freq(r._has_minim_freq),
  _magic(r._magic),
  _freq_order(nullptr),
  _repart_table(r._repart_table),
  is_load(r.is_load)
{
  if (r._freq_order)
  {
    _freq_order = new uint32_t[_nb_minims];
    copy(&r._freq_order[0], &r._freq_order[_nb_minims-1], _freq_order);
  }
}


RepartFile& RepartFile::operator=(const RepartFile &r)
{
  _path = r._path;
  _path_freq = r._path_freq;
  _nb_part = r._nb_part;
  _nb_minims = r._nb_minims;
  _has_minim_freq = r._has_minim_freq;
  _magic = r._magic;
  _freq_order = nullptr;
  _repart_table = r._repart_table;
  is_load = r.is_load;

  if (r._freq_order)
  {
    _freq_order = new uint32_t[_nb_minims];
    copy(&r._freq_order[0], &r._freq_order[_nb_minims-1], _freq_order);
  }
  return *this;
}


void RepartFile::load()
{
  ifstream fin(_path, ios::binary | ios::in);
  fin.read((char*)&_nb_part, sizeof(_nb_part));
  fin.read((char*)&_nb_minims, sizeof(_nb_minims));
  fin.read((char*)&_nb_pass, sizeof(_nb_pass));
  _repart_table.resize(_nb_minims);

  fin.read((char*)_repart_table.data(), sizeof(uint16_t)*_nb_minims);

  fin.read((char*)&_has_minim_freq, sizeof(bool));
  fin.read((char*)&_magic, sizeof(_magic));
  if (_magic != MAGIC_NUMBER)
    throw runtime_error("Unable to load " + _path + ", possibly due to bad format.");
  fin.close();

  if (_has_minim_freq && !_path_freq.empty())
  {
    ifstream fin(_path_freq, ios::binary | ios::in);
    _freq_order = new uint32_t[_nb_minims];
    fin.read((char*)_freq_order, sizeof(uint32_t)*_nb_minims);
    fin.read((char*)&_magic, sizeof(_magic));
    if (_magic != MAGIC_NUMBER)
      throw runtime_error("Unable to load " + _path + ", possibly due to bad format.");
    fin.close();
  }
  else
  {
    _freq_order = nullptr;
  }
  is_load = true;
}


uint16_t RepartFile::get(uint64_t minimizer_value)
{
  if (!is_load)
    load();
  return _repart_table[minimizer_value];
}
#endif
} // end of namespace km