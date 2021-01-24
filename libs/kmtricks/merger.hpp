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
#include <stdexcept>
#include <exception>
#include <unordered_map>

#include "utilities.hpp"
#include "sequences.hpp"
#include "lz4_stream.hpp"
#include "io.hpp"

typedef unsigned int  uint;
typedef unsigned char uchar;

using namespace std;

namespace km
{

template<typename K, typename C, typename F>
class Merger
{
public:

  struct hshcount_t
  {
    K khash;
    C count;
    bool khash_set;
  };

  struct stream_t
  {
    int eof, iskhash;
    F *f;
  };

public:
  Merger(string &fof_path, uint abundance, uint recurrence, uint header_size, bool vector,  uint save_if=0, bool stats=false);
  Merger(string &fof_path, vector<uint>& abundances, uint recurrence, uint header_size, bool vector, uint save_if=0, bool stats=false);

  ~Merger();

  Merger(const Merger<K, C, F> &m);

  Merger<K, C, F> &operator=(const Merger<K, C, F> &m);

  int init();

  void next();

  km::Kmer<K> get_kmer(size_t ksize);

private:
  inline void initp();

  inline int readb(size_t i);

public:
  vector<string>   pfiles;
  bool             keep;
  bool             end;
  K                m_khash;
  vector<C>        counts;
  size_t           nb_files;
  size_t           vlen;
  uchar            *_bit_vector;
  vector<uint64_t> _non_solid;
  vector<uint64_t> _saved;
  vector<uint>     total;
  vector<uint>     total_w_saved;

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

  vector<stream_t *>          _st;
  vector<hshcount_t *>  _hc;
  vector<uchar *>             _headers;
  vector<uint>                _abs_vec;
  vector<uint64_t>            _need_check;

  uint _save_if;
  bool _stats;
};

#ifndef _KM_LIB_INCLUDE_
template<typename K, typename C, typename F>
Merger<K, C, F>::Merger(string &fof_path, uint abundance, uint recurrence, uint header_size, bool vector, uint save_if, bool stats)
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
    _vector(vector),
    _save_if(save_if),
    _stats(stats)
{
  init();
}

template<typename K, typename C, typename F>
Merger<K, C, F>::Merger(string& fof_path, vector<uint>& abundances, uint recurrence, uint header_size, bool vector, uint save_if, bool stats)
  : keep(false),
    end(false),
    m_khash(0),
    nb_files(0),
    vlen(0),
    _bit_vector(nullptr),
    _path(fof_path),
    _a_min(0),
    _r_min(recurrence),
    _nm_khash(0),
    _buf_size((sizeof(K) + sizeof(C)) * 128),
    _hsize(header_size),
    _m_k_set(false),
    _nm_kh_set(false),
    _vector(vector),
    _abs_vec(abundances),
    _save_if(save_if),
    _stats(stats)
{
  init();
}


template<typename K, typename C, typename F>
Merger<K, C, F>::~Merger()
{
  for ( size_t i = 0; i < nb_files; i++ )
  {
    if ( _headers.size())
      delete[] _headers[i];
    if ( _st[i]->f ) delete _st[i]->f;
    delete _st[i];
    delete _hc[i];
  }
  if ( _bit_vector )
    delete[] _bit_vector;
  // useless
}


template<typename K, typename C, typename F>
Merger<K, C, F>::Merger(const Merger<K, C, F> &m)
  : keep(m.keep),
    end(m.end),
    m_khash(m.m_khash),
    nb_files(m.nb_files),
    vlen(m.vlen),
    _bit_vector(nullptr),
    //_bit_vector(m._bit_vector)
    _path(m._path),
    _a_min(m._a_min),
    _r_min(m._r_min),
    _nm_khash(m._nm_khash),
    _buf_size(m._buf_size),
    _hsize(m._hsize),
    _m_k_set(m._m_k_set),
    _nm_kh_set(m._nm_kh_set),
    _vector(m._vector),
    _abs_vec(m._abs_vec),
    _save_if(m._save_if),
    _stats(m._stats),
    _saved(m._saved),
    _non_solid(m._non_solid),
    _need_check(m._need_check),
    total(m.total),
    total_w_saved(m.total_w_saved)
{
  _hc.resize(nb_files);
  _st.resize(nb_files);

  for ( size_t i = 0; i < nb_files; i++ )
  {
    _hc[i] = new hshcount_t();
    _hc[i]->khash = m._hc[i]->khash;
    _hc[i]->count = m._hc[i]->count;
    _hc[i]->khash_set = m._hc[i]->khash_set;

    _st[i] = new stream_t();
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


template<typename K, typename C, typename F>
Merger<K, C, F> &Merger<K, C, F>::operator=(const Merger<K, C, F> &m)
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
  _abs_vec = m._abs_vec;
  _save_if = m._save_if;
  _stats = m._stats;
  _saved = m._saved;
  _non_solid = m._non_solid;
  _need_check = m._need_check;
  total = m.total;
  total_w_saved = m.total_w_saved;

  _hc.resize(nb_files);
  _st.resize(nb_files);

  for ( size_t i = 0; i < nb_files; i++ )
  {
    _hc[i] = new hshcount_t();
    _hc[i]->khash = m._hc[i]->khash;
    _hc[i]->count = m._hc[i]->count;
    _hc[i]->khash_set = m._hc[i]->khash_set;

    _st[i] = new stream_t();
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


template<typename K, typename C, typename F>
void Merger<K, C, F>::initp()
{
  ifstream fof(_path, ios::in);
  if (!fof.good())
    throw runtime_error("Unable to open " + _path + ".");
  string line;
  while ( getline(fof, line))
  {
    pfiles.push_back(line);
  }
  nb_files = pfiles.size();
}


template<typename K, typename C, typename F>
int Merger<K, C, F>::readb(size_t i)
{
  bool ret = _st[i]->f->read(_hc[i]->khash, _hc[i]->count);
  if (!ret)
  {
    _st[i]->eof = 1;
    return 0;
  }
  return 1;
}


template<typename K, typename C, typename F>
int Merger<K, C, F>::init()
{
  initp();
  if ( _vector )
  {
    vlen = NBYTE(nb_files);
    _bit_vector = new uchar[vlen]();
  }
  _m_k_set = false;
  end = false;

  counts.resize(nb_files, 0);
  total.resize(nb_files, 0);
  total_w_saved.resize(nb_files, 0);

  if (!_a_min)
    if(_abs_vec.size() > 0 && _abs_vec.size() != nb_files)
        throw runtime_error("Nb files != abundance_vec.size()");

  if (_save_if)
    _need_check.resize(nb_files);
  
  if (_stats)
  {
    _non_solid.resize(nb_files);
    _saved.resize(nb_files);
  }
  
  for ( size_t i = 0; i < nb_files; i++ )
  {
    _hc.push_back(new hshcount_t());
    _st.push_back(new stream_t());

    _st[i]->f = new F(pfiles[i]);

    //if ( _hsize > 0 )
    //{
    //  _headers.push_back(new uchar[_hsize]());
    //  _st[i]->f->read((char*)_headers[i], _hsize-1);
    //  _headers[i][_hsize - 1] = '\0';
    //}

    if ( !readb(i))
    {
      _hc[i]->khash_set = false;
      delete _st[i]->f;
      _st[i]->f = nullptr;
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


template<typename K, typename C, typename F>
void Merger<K, C, F>::next()
{
  uint rec = 0;
  uint solid_in = 0;
  keep = false;
  end = true;
  m_khash = _nm_khash;
  _nm_kh_set = false;
  _need_check.clear();
  if ( _bit_vector )
    memset(_bit_vector, 0, vlen);
  
  for ( size_t i = 0; i < nb_files; i++ )
  {
    if ( _hc[i]->khash_set && _hc[i]->khash == m_khash )
    {
      end = false;
      counts[i] = _hc[i]->count;

      if ( _a_min ? counts[i] >= _a_min : counts[i] >= _abs_vec[i] )
      {
        rec++;
        solid_in++;
        if ( _vector )
          BITSET(_bit_vector, i);
        total[i] += counts[i];
      }
      else
      {
        if (_stats)
          _non_solid[i]++;
        if (_save_if)
          _need_check.push_back(i);
        else
          counts[i] = 0;
      }

      if ( !readb(i))
      {
        _hc[i]->khash_set = false;
        delete _st[i]->f;
        _st[i]->f = nullptr;
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
  for (auto& p : _need_check)
  {
    if (solid_in >= _save_if)
    {
      _saved[p]++;
      total_w_saved[p] += counts[p];
      if (_vector)
        BITSET(_bit_vector, p);
    }
    else
      counts[p] = 0;
  }
  if ( rec >= _r_min )
    keep = true;
}


template<typename K, typename C, typename F>
km::Kmer<K> Merger<K, C, F>::get_kmer(size_t ksize)
{
  return km::Kmer<K>(m_khash, ksize, false);
}
#endif
} // end of namespace km