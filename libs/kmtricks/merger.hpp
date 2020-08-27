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
#include <cstring>
#include <iostream>
#include <exception>
#include <zlib.h>

#include "utilities.hpp"
#include "sequences.hpp"

typedef unsigned int  uint;
typedef unsigned char uchar;

using namespace std;

namespace km
{

struct stream_t
{
  uchar *buf;
  int begin, end, eof, iskhash;
  gzFile f;
};


template<typename K, typename C>
struct hshcount_t
{
  K khash;
  C count;
  bool khash_set;
};


template<typename K, typename C>
class Merger
{
public:
  Merger(string &fof_path, uint abundance, uint recurrence, uint header_size, bool vector);

  ~Merger();

  Merger(const Merger<K, C> &m);

  Merger<K, C> &operator=(const Merger<K, C> &m);

  int init();

  void next();

  km::Kmer<K> get_kmer(size_t ksize);

private:
  inline void initp();

  inline int readb(size_t i);

public:
  vector<string> pfiles;
  bool keep;
  bool end;
  K m_khash;
  vector<C> counts;
  size_t nb_files;
  size_t vlen;
  uchar *_bit_vector;

private:
  string _path;
  uint _a_min;
  uint _r_min;

  K _nm_khash;

  size_t _buf_size;
  size_t _hsize;

  bool _m_k_set;
  bool _nm_kh_set;
  bool _vector;

  vector<stream_t *> _st;
  vector<hshcount_t<K, C> *> _hc;
  vector<uchar *> _headers;
};

#ifndef _KM_LIB_INCLUDE_
template<typename K, typename C>
Merger<K, C>::Merger(string &fof_path, uint abundance, uint recurrence, uint header_size, bool vector)
  : keep(false),
    end(false),
    m_khash(0),
    nb_files(0),
    vlen(0),
    _bit_vector(nullptr),
    _path(fof_path),
    _a_min(abundance),
    _r_min(recurrence),
    _nm_khash(0),
    _buf_size((sizeof(K) + sizeof(C)) * 128),
    _hsize(header_size),
    _m_k_set(false),
    _nm_kh_set(false),
    _vector(vector)
{
  init();
}


template<typename K, typename C>
Merger<K, C>::~Merger()
{
  for ( size_t i = 0; i < nb_files; i++ )
  {
    if ( _headers.size())
      delete[] _headers[i];
    delete _st[i];
    delete _hc[i];
  }
  if ( _bit_vector )
    delete[] _bit_vector;
}


template<typename K, typename C>
Merger<K, C>::Merger(const Merger<K, C> &m)
  : keep(m.keep),
    end(m.end),
    m_khash(m.m_khash),
    nb_files(m.nb_files),
    vlen(m.vlen),
    _bit_vector(nullptr),
    _path(m._path),
    _a_min(m._a_min),
    _r_min(m._r_min),
    _nm_khash(m._nm_khash),
    _buf_size(m._buf_size),
    _hsize(m._hsize),
    _m_k_set(m._m_k_set),
    _nm_kh_set(m._nm_kh_set),
    _vector(m._vector)
{
  _hc.resize(nb_files);
  _st.resize(nb_files);

  for ( size_t i = 0; i < nb_files; i++ )
  {
    _hc[i] = new hshcount_t<K, C>();
    _hc[i]->khash = m._hc[i]->khash;
    _hc[i]->count = m._hc[i]->count;
    _hc[i]->khash_set = m._hc[i]->khash_set;

    _st[i] = new stream_t();
    _st[i]->buf = m._st[i]->buf;
    _st[i]->begin = m._st[i]->begin;
    _st[i]->end = m._st[i]->end;
    _st[i]->eof = m._st[i]->eof;
    _st[i]->iskhash = m._st[i]->iskhash;
    _st[i]->f = m._st[i]->f;
  }

  if ( m._headers.size())
  {
    _headers.resize(nb_files);
    for ( size_t i = 0; i < nb_files; i++ )
    {
      _headers[i] = new uchar[_hsize]();
      copy(&(_headers[i][0]), &(_headers[i][_hsize - 1]), _headers[i]);
    }
  }

  if ( m._bit_vector )
  {
    _bit_vector = new uchar[vlen]();
    copy(&m._bit_vector[0], &m._bit_vector[vlen - 1], _bit_vector);
  }
}


template<typename K, typename C>
Merger<K, C> &Merger<K, C>::operator=(const Merger<K, C> &m)
{
  keep = m.keep;
  end = m.end;
  m_khash = m.m_khash;
  nb_files = m.nb_files;
  vlen = m.vlen;
  _path = m._path;
  _a_min = m._a_min;
  _r_min = m._r_min;
  _nm_khash = m._nm_khash;
  _buf_size = m._buf_size;
  _hsize = m._hsize;
  _m_k_set = m._m_k_set;
  _nm_kh_set = m._nm_kh_set;
  _vector = m._vector;

  _hc.resize(nb_files);
  _st.resize(nb_files);

  for ( size_t i = 0; i < nb_files; i++ )
  {
    _hc[i] = new hshcount_t<K, C>();
    _hc[i]->khash = m._hc[i]->khash;
    _hc[i]->count = m._hc[i]->count;
    _hc[i]->khash_set = m._hc[i]->khash_set;

    _st[i] = new stream_t();
    _st[i]->buf = m._st[i]->buf;
    _st[i]->begin = m._st[i]->begin;
    _st[i]->end = m._st[i]->end;
    _st[i]->eof = m._st[i]->eof;
    _st[i]->iskhash = m._st[i]->iskhash;
    _st[i]->f = m._st[i]->f;
  }

  if ( m._headers.size())
  {
    _headers.resize(nb_files);
    for ( size_t i = 0; i < nb_files; i++ )
    {
      _headers[i] = new uchar[_hsize]();
      copy(&(_headers[i][0]), &(_headers[i][_hsize - 1]), _headers[i]);
    }
  }

  if ( m._bit_vector )
  {
    _bit_vector = new uchar[vlen]();
    copy(&m._bit_vector[0], &m._bit_vector[vlen - 1], _bit_vector);
  }
  return *this;
}


template<typename K, typename C>
void Merger<K, C>::initp()
{
  ifstream fof(_path, ios::in);
  string line;
  while ( getline(fof, line))
  {
    pfiles.push_back(line);
  }
  nb_files = pfiles.size();
}


template<typename K, typename C>
int Merger<K, C>::readb(size_t i)
{
  if ( _st[i]->begin >= _st[i]->end )
  {
    _st[i]->begin = 0;
    _st[i]->end = gzread(_st[i]->f, _st[i]->buf, _buf_size);
    if ( _st[i]->end == 0 )
    {
      _st[i]->eof = 1;
      return 0;
    }
  }
  memcpy(&_hc[i]->khash, &_st[i]->buf[_st[i]->begin], sizeof(K));
  _st[i]->begin += sizeof(K);

  memcpy(&_hc[i]->count, &_st[i]->buf[_st[i]->begin], sizeof(C));
  _st[i]->begin += sizeof(C);

  return 1;
}


template<typename K, typename C>
int Merger<K, C>::init()
{
  initp();
  if ( _vector )
  {
    vlen = NMOD8(NBYTE(nb_files));
    _bit_vector = new uchar[vlen]();
  }
  _m_k_set = false;
  end = false;
  counts.reserve(nb_files);

  for ( size_t i = 0; i < nb_files; i++ )
  {
    _hc.push_back(new hshcount_t<K, C>());
    _st.push_back(new stream_t());

    _st[i]->f = gzopen(pfiles[i].c_str(), "r");
    if ( !_st[i]->f )
      throw runtime_error("Unable to open " + pfiles[i]);

    _st[i]->buf = new uchar[_buf_size]();

    if ( _hsize > 0 )
    {
      _headers.push_back(new uchar[_hsize]());
      gzread(_st[i]->f, _headers[i], _hsize - 1);
      _headers[i][_hsize - 1] = '\0';
    }

    if ( !readb(i))
    {
      _hc[i]->khash_set = false;
      gzclose(_st[i]->f);
      delete[] _st[i]->buf;
    }
    else
      _hc[i]->khash_set = true;

    if (((!_m_k_set) || _hc[i]->khash < m_khash) && _hc[i]->khash_set )
    {
      m_khash = _nm_khash = _hc[i]->khash;
      _m_k_set = true;
    }
  }
  return 1;
}


template<typename K, typename C>
void Merger<K, C>::next()
{
  uint rec = 0;
  keep = false;
  end = true;
  m_khash = _nm_khash;
  _nm_kh_set = false;

  if ( _bit_vector )
    memset(_bit_vector, 0, vlen);
  for ( size_t i = 0; i < nb_files; i++ )
  {
    if ( _hc[i]->khash_set && _hc[i]->khash == m_khash )
    {
      end = false;
      counts[i] = _hc[i]->count;

      if ( _vector )
        BITSET(_bit_vector, i);

      if ( counts[i] >= _a_min )
        rec++;

      if ( !readb(i))
      {
        _hc[i]->khash_set = false;
        gzclose(_st[i]->f);
        delete[] _st[i]->buf;
      }
    }
    else
      counts[i] = 0;

    if ( _hc[i]->khash_set && (!_nm_kh_set || _hc[i]->khash < _nm_khash))
    {
      _nm_khash = _hc[i]->khash;
      _nm_kh_set = true;
    }
  }
  if ( rec >= _r_min )
    keep = true;
}


template<typename K, typename C>
km::Kmer<K> Merger<K, C>::get_kmer(size_t ksize)
{
  return km::Kmer<K>(m_khash, ksize, false);
}
#endif
} // end of namespace km