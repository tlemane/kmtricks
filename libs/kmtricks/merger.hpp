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

typedef unsigned int  uint;
typedef unsigned char uchar;

#define NBYTE(bits) (((bits) >> 3) + ((bits) % 8 != 0))
#define NMOD8(byte) ((byte)+(8-((byte)%8)))

#define BITMASK(b) (1 << ((b) % CHAR_BIT))
#define BITSLOT(b) ((b) / CHAR_BIT)
#define BITSET(a, b) ((a)[BITSLOT(b)] |= BITMASK(b))

using namespace std;

static const char bToN[] = {'A', 'C', 'T', 'G'};

struct stream_t
{
  uchar* buf;
  int begin, end, eof, iskhash;
  gzFile f;
};

template<unsigned long long klength>
struct requiredK
{
  enum {value =
    klength <= 4  ?  8 :
    klength <= 8  ? 16 :
    klength <= 16 ? 32 :
    klength <= 32 ? 64 :
                    128
  };
};

template<unsigned long long klength>
struct requiredC
{
  enum {value =
    klength <= 0xFF       ?  8 :
    klength <= 0xFFFF     ? 16 :
                            32
  };
};
template<int bits> struct select_;
template<> struct select_ <8> {typedef uint8_t type;};
template<> struct select_ <16> {typedef uint16_t type;};
template<> struct select_ <32> {typedef uint32_t type;};
template<> struct select_ <64> {typedef uint64_t type;};
template<> struct select_ <128> {typedef __uint128_t type;};

template<unsigned long long klength>
struct selectK : select_<requiredK<klength>::value> {};

template<unsigned long long klength>
struct selectC : select_<requiredC<klength>::value> {};

template <typename K, typename C>
struct hshcount_t
{
  K khash;
  C count;
  bool khash_set;
};

template <typename K, typename C>
class Merger
{
public:
  Merger(string &fof_path, uint abundance, uint recurrence, uint header_size, bool vector);
  int         init();
  void        next();
  string      getk(size_t ksize, const char* binToNt);
  void        destroy();

private:
  inline void initp();
  inline int  readb(size_t i);

public:
  vector<string>  pfiles;
  bool            keep;
  bool            end;
  K               m_khash;
  vector<C>       counts;
  size_t          nb_files;
  size_t          vlen;
  uchar           *_bit_vector;

private:
  uint _a_min;
  uint _r_min;
  string _path;

  K _nm_khash;

  size_t _buf_size;
  size_t _hsize;

  bool _m_k_set;
  bool _nm_kh_set;
  bool _vector;

  vector<stream_t*>         _st;
  vector<hshcount_t<K, C>*> _hc;
  vector<uchar*>            _headers;

};

template <typename K, typename C>
Merger<K, C>::Merger(string &fof_path, uint abundance, uint recurrence, uint header_size, bool vector)
  : _path(fof_path), _a_min(abundance), _r_min(recurrence), _hsize(header_size), _vector(vector)
{
  _buf_size = (sizeof(K)+sizeof(C))*128;
  init();
}

template <typename K, typename C>
void Merger<K, C>::initp()
{
  ifstream fof(_path, ios::in);
  string line;
  while (getline(fof, line))
  {
    pfiles.push_back(line);
  }
  nb_files = pfiles.size();
}

template <typename K, typename C>
int Merger<K, C>::readb(size_t i)
{
  if (_st[i]->begin >= _st[i]->end)
  {
    _st[i]->begin = 0;
    _st[i]->end = gzread(_st[i]->f, _st[i]->buf, _buf_size);
    if (_st[i]->end == 0) { _st[i]->eof = 1; return 0; }
  }
  memcpy(&_hc[i]->khash, &_st[i]->buf[_st[i]->begin], sizeof(K));
  _st[i]->begin += sizeof(K);

  memcpy(&_hc[i]->count, &_st[i]->buf[_st[i]->begin], sizeof(C));
  _st[i]->begin += sizeof(C);

  return 1;
}

template <typename K, typename C>
int Merger<K, C>::init()
{
  initp();
  if (_vector)
  {
    vlen = NMOD8(NBYTE(nb_files));
    _bit_vector = new uchar[vlen]();
  }
  _m_k_set  = false;
  end       = false;
  counts.reserve(nb_files);

  for (size_t i=0; i<nb_files; i++)
  {
    _hc.push_back(new hshcount_t<K, C>());
    _st.push_back(new stream_t());

    _st[i]->f = gzopen(pfiles[i].c_str(), "r");
    if (!_st[i]->f)
      throw runtime_error("Unable to open " + pfiles[i]);

    _st[i]->buf = new uchar[_buf_size]();

    if (_hsize > 0)
    {
      _headers.push_back(new uchar[_hsize]());
      gzread(_st[i]->f, _headers[i], _hsize-1);
      _headers[i][_hsize-1] = '\0';
    }

    if ( !readb(i) )
    {
      _hc[i]->khash_set = false;
      gzclose(_st[i]->f);
      delete[] _st[i]->buf;
    }
    else
      _hc[i]->khash_set = true;

    if ( ( (!_m_k_set) || _hc[i]->khash < m_khash) && _hc[i]->khash_set)
    {
      m_khash = _nm_khash = _hc[i]->khash;
      _m_k_set = true;
    }
  }
  return 1;
}

template <typename K, typename C>
void Merger<K,C>::next()
{
  uint rec = 0;
  keep = false;
  end = true;
  m_khash = _nm_khash;
  _nm_kh_set = false;

  if (_bit_vector) memset(_bit_vector, 0, vlen);
  for (size_t i=0; i<nb_files; i++ )
  {
    if (_hc[i]->khash_set && _hc[i]->khash == m_khash)
    {
      end = false;
      counts[i] = _hc[i]->count;

      if (_vector) BITSET(_bit_vector, i);

      if (counts[i] >= _a_min) rec++;

      if ( !readb(i) )
      {
        _hc[i]->khash_set = false;
        gzclose(_st[i]->f);
        delete[] _st[i]->buf;
      }
    }
    else counts[i] = 0;

    if (_hc[i]->khash_set && (!_nm_kh_set || _hc[i]->khash < _nm_khash))
    {
      _nm_khash = _hc[i]->khash;
      _nm_kh_set = true;
    }
  }
  if (rec >= _r_min) keep = true;
}

template <typename K, typename C>
string Merger<K, C>::getk(size_t ksize, const char* binToNt)
{
  char tmp[ksize+1];
  K value = m_khash;
  for (int i=ksize-1; i>=0; i--)
  {
    tmp[i] = bToN[value&3];
    value = value >> 2;
  }
  tmp[ksize] = '\0';
  return tmp;
}

template<typename K, typename C>
void Merger<K, C>::destroy()
{
  for (size_t i=0; i<nb_files; i++)
  {
    if (_headers.size())
      delete[] _headers[i];
    delete _st[i];
    delete _hc[i];
  }
  if (_bit_vector) delete[] _bit_vector;
}

