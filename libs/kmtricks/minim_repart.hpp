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
static const uint32_t MAGIC_NUMBER = 0x12345678;
static const char bToN[] = {'A', 'C', 'T', 'G'};
static const char rev[] = {'T', 'G', 'A', 'C'};

static const uint8_t NToB[256] = {
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 0, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 0, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};



// from https://github.com/tlemane/gatb-core-stripped/blob/7ecaff35b97ce2bfb7731dc2e6dbc49719e23a98/src/gatb/kmer/impl/Model.hpp#L1220
bool is_allowed (uint64_t mmer, uint32_t len)
		{	
			u_int64_t  _mmask_m1  ;
			u_int64_t  _mask_0101 ;
			u_int64_t  _mask_ma1 ;
			
			//code to ban mmer with AA inside except if at the beginnning
			// A C T G        00   01   10   11
			_mmask_m1  = (1 << ((len-2)*2)) -1 ; //removes two first letters m = 8  gives    00 00 11 11 11 11 11 11
			_mask_0101 = 0x5555555555555555  ; //         01 01 01 01 01 01 01 01
			_mask_ma1  = _mask_0101 & _mmask_m1;//        00 00 01 01 01 01 01 01
			
			u_int64_t a1 = mmer; //
			a1 =   ~(( a1 )   | (  a1 >>2 ));  //
			a1 =((a1 >>1) & a1) & _mask_ma1 ;  //
			
			if(a1 != 0) return false;
			
			// if ((mmer & 0x3f) == 0x2a)   return false;   // TTT suffix
			// if ((mmer & 0x3f) == 0x2e)   return false;   // TGT suffix
			// if ((mmer & 0x3c) == 0x28)   return false;   // TT* suffix
			// for (uint32_t j = 0; j < len - 3; ++j)       // AA inside
			//      if ((mmer & 0xf) == 0)  return false;
			//      else                    mmer >>= 2;
			// if (mmer == 0)               return false;   // AAA prefix
			// if (mmer == 0x04)            return false;   // ACA prefix
			// if ((mmer & 0xf) == 0)   return false;       // *AA prefix
			
			return true;
		}


class RepartFile
{
public:
  RepartFile(string m_path, string f_path = "");
  void load();
  uint16_t get(uint64_t minim_value);

public:
  bool _is_load;
private:
  string   _path;
  string   _path_freq;
  uint16_t _nb_part;
  uint64_t _nb_minims;
  uint16_t _nb_pass;
  bool     _has_minim_freq;
  uint32_t _magic;
  vector<uint16_t> _repart_table;
  uint32_t* _freq_order;
};

RepartFile::RepartFile(string m_path, string f_path)
  : _path(m_path), _path_freq(f_path), _is_load(false)
{ }

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

  if(_has_minim_freq)
  {
    ifstream fin(_path_freq, ios::binary | ios::in);
    _freq_order = new uint32_t[_nb_minims];
    fin.read((char*)_freq_order, sizeof(uint32_t)*_nb_minims);
    fin.read((char*)&_magic, sizeof(_magic));
    if (_magic != MAGIC_NUMBER)
      throw runtime_error("Unable to load " + _path + ", possibly due to bad format.");
  }
  else
  {
    _freq_order = nullptr;
  }
  _is_load = true;
}

uint16_t RepartFile::get(uint64_t minim_value)
{
  return _repart_table[minim_value];
}

template<typename KT>
class MinimRepart
{
public:
  MinimRepart(RepartFile& rfile);
  uint64_t  get_minim_from_str(string seq, size_t s_size, size_t m_size);
  string    int_to_str(KT seq, size_t size);
  KT        seq_to_int(string seq, size_t s_size);
  KT        rev_comp(KT seq, size_t size);

  uint16_t get_partition(uint64_t minim_value);


private:
  RepartFile  _rfile;
};

template<typename KT>
MinimRepart<KT>::MinimRepart(RepartFile &rfile)
  : _rfile(rfile)
{
  if (!_rfile._is_load) _rfile.load();
}

template<typename KT>
uint16_t MinimRepart<KT>::get_partition(uint64_t minim_value)
{
  return _rfile.get(minim_value);
}

template<typename KT>
KT MinimRepart<KT>::seq_to_int(string seq, size_t s_size)
{
  KT res = 0;
  for (size_t c=0; c<s_size; c++)
  {
    res <<= 2;
    res |= NToB[seq[c]];
  }
  return res;
}

template<typename KT>
uint64_t MinimRepart<KT>::get_minim_from_str(string seq, size_t s_size, size_t m_size)
{
  KT kmer = seq_to_int(seq, s_size);
  KT revcomp = rev_comp(kmer, s_size);
  kmer = revcomp < kmer ? revcomp : kmer;
  uint64_t minim = numeric_limits<uint64_t>::max();
  uint64_t _mmer_mask = numeric_limits<uint64_t>::max() >> (64-(m_size*2));
  int nb_mmers = s_size - m_size + 1;
  uint64_t tmp;
  for (int i=nb_mmers-1; i>=0; i--)
  {
    tmp = (uint64_t)(kmer >> (i*2)) & _mmer_mask;

    if (is_allowed(tmp, m_size) && tmp < minim) minim = tmp;
  }
  return minim;
}

template<typename KT>
string MinimRepart<KT>::int_to_str(KT seq, size_t size)
{
  char tmp[size+1];
  KT value = seq;
  for (int i=size-1; i>=0; i--)
  {
    tmp[i] = bToN[value&3];
    value = value >> 2;
  }
  tmp[size] = '\0';
  return tmp;
}

template<typename KT>
KT MinimRepart<KT>::rev_comp(KT seq, size_t size)
{
  KT res = 0;
  for (int c=size-1; c>=0; c--)
  {
    res <<= 2;
    res |= NToB[rev[seq&3]];
    seq >>= 2;
  }
  KT _kmer_mask = numeric_limits<KT>::max() >> ((sizeof(KT)*8)-(size*2));
  return res & _kmer_mask;
}

