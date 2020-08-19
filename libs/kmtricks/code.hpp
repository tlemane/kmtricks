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
#include <cstring>

using namespace std;

typedef unsigned char uchar;

uchar bToN[] = {'A', 'C', 'T', 'G'};
uchar revC[] = {'T', 'G', 'A', 'C'};

uchar NToB[256] = {
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

namespace km
{

template<typename K>
class Code
{
public:
  Code();

  Code(uchar *bits_to_nt, uchar *nt_to_bits, uchar *revc);

  Code(uchar *bits_to_nt);

  ~Code();

  Code(const Code<K> &c);

  Code<K> &operator=(const Code<K> &c);

  void set_default_encoding();

  void set_encoding(uchar *bits_to_nt, uchar *nt_to_bits, uchar *revc);

  void set_encoding(uchar *bits_to_nt);

  K encode(string value, size_t size);

  uchar encode(char value);

  string decode(K value, size_t size);

  string decode(uchar value);

private:
  void destroy();

public:
  uchar *_bToN;
  uchar *_NToB;
  uchar *_revC;

private:
  bool _custom_enc;
  bool _from_out;
};

template<typename K>
Code<K>::~Code()
{
  if ( _custom_enc && !_from_out )
    destroy();
}

template<typename K>
void Code<K>::destroy()
{
  delete[] _bToN;
  delete[] _NToB;
  delete[] _revC;
}

template<typename K>
Code<K>::Code()
  :
  _bToN(nullptr),
  _NToB(nullptr),
  _revC(nullptr),
  _custom_enc(false),
  _from_out(false)
{
  set_default_encoding();
}

template<typename K>
Code<K>::Code(uchar *bits_to_nt, uchar *nt_to_bits, uchar *revc)
  :
  _bToN(nullptr),
  _NToB(nullptr),
  _revC(nullptr),
  _custom_enc(true),
  _from_out(true)
{
  set_encoding(bits_to_nt, nt_to_bits, revc);
}

template<typename K>
Code<K>::Code(uchar *bits_to_nt)
  :
  _bToN(nullptr),
  _NToB(nullptr),
  _revC(nullptr),
  _custom_enc(true),
  _from_out(false)
{
  set_encoding(bits_to_nt);
}

template<typename K>
Code<K>::Code(const Code<K> &c)
  :
  _bToN(nullptr),
  _NToB(nullptr),
  _revC(nullptr),
  _custom_enc(c._custom_enc),
  _from_out(c._from_out)
{
  if ( _custom_enc && !_from_out )
  {
    _bToN = new uchar[4]();
    _NToB = new uchar[256]();
    _revC = new uchar[4]();

    copy(&c._bToN[0], &c._bToN[3], _bToN);
    copy(&c._revC[0], &c._revC[3], _revC);
    copy(&c._NToB[0], &c._NToB[255], _NToB);
  }
  else
  {
    _bToN = c._bToN;
    _NToB = c._NToB;
    _revC = c._revC;
  }
}

template<typename K>
Code<K> &Code<K>::operator=(const Code<K> &c)
{
  _custom_enc = c._custom_enc;
  _from_out = c._from_out;
  if ( _custom_enc && !_from_out )
  {
    _bToN = new uchar[4]();
    _NToB = new uchar[256]();
    _revC = new uchar[4]();

    copy(&c._bToN[0], &c._bToN[3], _bToN);
    copy(&c._revC[0], &c._revC[3], _revC);
    copy(&c._NToB[0], &c._NToB[255], _NToB);
  }
  else
  {
    _bToN = c._bToN;
    _NToB = c._NToB;
    _revC = c._revC;
  }
  return *this;
}

template<typename K>
void Code<K>::set_default_encoding()
{
  _bToN = bToN;
  _NToB = NToB;
  _revC = revC;
  _custom_enc = false;
}

template<typename K>
void Code<K>::set_encoding(uchar *bits_to_nt, uchar *nt_to_bits, uchar *revc)
{
  if ( !_custom_enc )
    destroy();
  _custom_enc = true;
  _from_out = true;
  _bToN = bits_to_nt;
  _NToB = nt_to_bits;
  _revC = revc;
}

template<typename K>
void Code<K>::set_encoding(uchar *bits_to_nt)
{
  if ( _custom_enc && !_from_out )
    destroy();
  _bToN = new uchar[4]();
  _NToB = new uchar[256]();
  _revC = new uchar[4]();

  char A, C, G, T;

  memcpy(_bToN, bits_to_nt, sizeof(uchar) * 4);
  for ( uchar c = 0; c < 4; c++ )
  {
    switch ( _bToN[c] )
    {
      case 'A':
        _revC[c] = 'T';
        A = c;
        break;
      case 'C':
        _revC[c] = 'G';
        C = c;
        break;
      case 'G':
        _revC[c] = 'C';
        G = c;
        break;
      case 'T':
        _revC[c] = 'A';
        T = c;
        break;
    }
  }

  for ( int i = 0; i < 256; i++ )
  {
    switch ( i )
    {
      case 'A':
      case 'a':
        _NToB[i] = A;
        break;
      case 'C':
      case 'c':
        _NToB[i] = C;
        break;
      case 'G':
      case 'g':
        _NToB[i] = G;
        break;
      case 'T':
      case 't':
        _NToB[i] = T;
        break;
    }
  }
  _custom_enc = true;
  _from_out = false;
}

template<typename K>
K Code<K>::encode(string value, size_t size)
{
  K res = 0;
  for ( size_t c = 0; c < size; c++ )
  {
    res <<= 2;
    res |= _NToB[(uchar) value[c]];
  }
  return res;
}

template<typename K>
uchar Code<K>::encode(char value)
{
  return _NToB[(uchar) value];
}

template<typename K>
string Code<K>::decode(K value, size_t size)
{
  char tmp[size + 1];
  for ( int i = size - 1; i >= 0; i-- )
  {
    tmp[i] = _bToN[value & 3];
    value = value >> 2;
  }
  tmp[size] = '\0';
  return tmp;
}

template<typename K>
string Code<K>::decode(uchar value)
{
  char tmp[5];
  for ( int i = 3; i >= 0; i-- )
  {
    tmp[i] = _bToN[value & 3];
    value = value >> 2;
  }
  tmp[4] = '\0';
  return tmp;
}

} // end of namespace