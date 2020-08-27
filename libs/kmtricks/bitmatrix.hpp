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
#include <exception>
#include <fstream>
#include <cstring>
#include <emmintrin.h>
#include <iostream>
#include <cassert>
#include <iomanip>
#include "utilities.hpp"

using namespace std;
typedef unsigned char uchar;

namespace km
{

#ifdef _KM_LIB_INCLUDE_
extern uchar reverseb[];
#endif

#ifndef _KM_LIB_INCLUDE_
uchar reverseb[] = {
  0x00, 0x80, 0x40, 0xc0, 0x20, 0xa0, 0x60, 0xe0,
  0x10, 0x90, 0x50, 0xd0, 0x30, 0xb0, 0x70, 0xf0,
  0x08, 0x88, 0x48, 0xc8, 0x28, 0xa8, 0x68, 0xe8,
  0x18, 0x98, 0x58, 0xd8, 0x38, 0xb8, 0x78, 0xf8,
  0x04, 0x84, 0x44, 0xc4, 0x24, 0xa4, 0x64, 0xe4,
  0x14, 0x94, 0x54, 0xd4, 0x34, 0xb4, 0x74, 0xf4,
  0x0c, 0x8c, 0x4c, 0xcc, 0x2c, 0xac, 0x6c, 0xec,
  0x1c, 0x9c, 0x5c, 0xdc, 0x3c, 0xbc, 0x7c, 0xfc,
  0x02, 0x82, 0x42, 0xc2, 0x22, 0xa2, 0x62, 0xe2,
  0x12, 0x92, 0x52, 0xd2, 0x32, 0xb2, 0x72, 0xf2,
  0x0a, 0x8a, 0x4a, 0xca, 0x2a, 0xaa, 0x6a, 0xea,
  0x1a, 0x9a, 0x5a, 0xda, 0x3a, 0xba, 0x7a, 0xfa,
  0x06, 0x86, 0x46, 0xc6, 0x26, 0xa6, 0x66, 0xe6,
  0x16, 0x96, 0x56, 0xd6, 0x36, 0xb6, 0x76, 0xf6,
  0x0e, 0x8e, 0x4e, 0xce, 0x2e, 0xae, 0x6e, 0xee,
  0x1e, 0x9e, 0x5e, 0xde, 0x3e, 0xbe, 0x7e, 0xfe,
  0x01, 0x81, 0x41, 0xc1, 0x21, 0xa1, 0x61, 0xe1,
  0x11, 0x91, 0x51, 0xd1, 0x31, 0xb1, 0x71, 0xf1,
  0x09, 0x89, 0x49, 0xc9, 0x29, 0xa9, 0x69, 0xe9,
  0x19, 0x99, 0x59, 0xd9, 0x39, 0xb9, 0x79, 0xf9,
  0x05, 0x85, 0x45, 0xc5, 0x25, 0xa5, 0x65, 0xe5,
  0x15, 0x95, 0x55, 0xd5, 0x35, 0xb5, 0x75, 0xf5,
  0x0d, 0x8d, 0x4d, 0xcd, 0x2d, 0xad, 0x6d, 0xed,
  0x1d, 0x9d, 0x5d, 0xdd, 0x3d, 0xbd, 0x7d, 0xfd,
  0x03, 0x83, 0x43, 0xc3, 0x23, 0xa3, 0x63, 0xe3,
  0x13, 0x93, 0x53, 0xd3, 0x33, 0xb3, 0x73, 0xf3,
  0x0b, 0x8b, 0x4b, 0xcb, 0x2b, 0xab, 0x6b, 0xeb,
  0x1b, 0x9b, 0x5b, 0xdb, 0x3b, 0xbb, 0x7b, 0xfb,
  0x07, 0x87, 0x47, 0xc7, 0x27, 0xa7, 0x67, 0xe7,
  0x17, 0x97, 0x57, 0xd7, 0x37, 0xb7, 0x77, 0xf7,
  0x0f, 0x8f, 0x4f, 0xcf, 0x2f, 0xaf, 0x6f, 0xef,
  0x1f, 0x9f, 0x5f, 0xdf, 0x3f, 0xbf, 0x7f, 0xff,
};
#endif

void __sse_trans(uint8_t const *inp, uint8_t *out, int nrows, int ncols);

class BitMatrix
{
public:
  BitMatrix(size_t n, size_t m, bool lendian, bool def);

  BitMatrix(string &file, size_t n, size_t m, bool lendian);

  BitMatrix(uchar *mat, size_t n, size_t m, bool lendian);

  ~BitMatrix();

  BitMatrix(const BitMatrix &b);

  BitMatrix &operator=(const BitMatrix &b);

  // i, j in bits
  void set_bit(size_t i, size_t j, bool value);

  void tog_bit(size_t i, size_t j);

  bool get_bit(size_t n, size_t m);

  // i, j in bytes
  void set_byte(size_t i, size_t j, uchar value);

  void tog_byte(size_t i, size_t j);

  uchar get_byte(size_t i, size_t j);

  uchar *get_line(size_t i);

  uchar *get_cols(size_t j);

  void clear();

  void dump(string file);

  void print_bytes();

  void print_bits();

  BitMatrix *transpose();

public:
  uchar *matrix;

private:
  inline void check8();

private:
  string _fpath;
  size_t _n;  // in bytes
  size_t _m;  // in bytes
  size_t _nb; // in bits
  size_t _mb; // in bits
  bool _le; // true if little endian
};

#ifndef _KM_LIB_INCLUDE_
void BitMatrix::check8()
{
  if ( _nb % 8 != 0 )
    throw runtime_error("n % 8 != 0 -> n=" + to_string(_nb));
  else if ( _mb % 8 != 0 )
    throw runtime_error("m % 8 != 0 -> m=" + to_string(_mb));
}


BitMatrix::~BitMatrix()
{
  if ( matrix )
    delete[] matrix;
}


BitMatrix::BitMatrix(size_t n, size_t m, bool lendian, bool def = false)
  : _n(n / 8),
    _m(m),
    _nb(n),
    _mb(m * 8),
    _le(lendian)
{
  check8();
  matrix = new uchar[_nb * _m]();
  if ( def )
    memset(matrix, 0xFF, _nb * _m);
}


BitMatrix::BitMatrix(string &file, size_t n, size_t m, bool lendian)
  : _fpath(file),
    _n(n / 8),
    _m(m),
    _nb(n),
    _mb(m * 8),
    _le(lendian)
{
  check8();
  matrix = new uchar[_nb * _m]();
  ifstream fin(_fpath, ios::in | ios::binary);
  if ( !fin )
    throw runtime_error("Unable to open : " + _fpath);
  fin.read((char *) matrix, _nb * _m);
  fin.close();
}


BitMatrix::BitMatrix(uchar *mat, size_t n, size_t m, bool lendian)
  : _n(n / 8),
    _m(m),
    _nb(n),
    _mb(m * 8),
    _le(lendian)
{
  check8();
  matrix = mat;
}


BitMatrix::BitMatrix(const BitMatrix &b)
  : _n(b._n),
    _m(b._m),
    _nb(b._nb),
    _mb(b._mb),
    _le(b._le)
{
  matrix = new uchar[_nb * _m]();
  copy(&b.matrix[0], &b.matrix[(_nb * _m) - 1], matrix);
}


BitMatrix &BitMatrix::operator=(const BitMatrix &b)
{
  _n = b._n;
  _m = b._m;
  _nb = b._nb;
  _mb = b._mb;
  _le = b._le;
  matrix = new uchar[_nb * _m]();
  copy(&b.matrix[0], &b.matrix[(_nb * _m) - 1], matrix);
  return *this;
}


void BitMatrix::set_bit(size_t i, size_t j, bool value)
{
  size_t offset = (i * _mb + j) / 8;
  size_t pos = (i * _mb + j) % 8;
  uchar mask = _le ? 0x1 << pos : 0x80 >> pos;
  if ( value )
    matrix[offset] |= mask;
  else
    matrix[offset] &= ~(mask);
}


void BitMatrix::tog_bit(size_t i, size_t j)
{
  size_t offset = (i * _mb + j) / 8;
  size_t pos = (i * _mb + j) % 8;
  uchar mask = _le ? 0x1 << pos : 0x80 >> pos;
  matrix[offset] ^= mask;
}


bool BitMatrix::get_bit(size_t i, size_t j)
{
  size_t offset = (i * _mb + j) / 8;
  size_t pos = (i * _mb + j) % 8;
  uchar mask = _le ? 0x1 << pos : 0x80 >> pos;
  return matrix[offset] & mask;
}


void BitMatrix::set_byte(size_t i, size_t j, uchar value) // i, j in bytes
{
  matrix[i * _m + j] = _le ? value : reverseb[value];
}


void BitMatrix::tog_byte(size_t i, size_t j)
{
  size_t offset = (i * _m + j);
  matrix[offset] = matrix[offset] ^ 0xFF;
}


uchar BitMatrix::get_byte(size_t i, size_t j)
{
  return matrix[i * _m + j];
}


uchar *BitMatrix::get_line(size_t i)
{
  uchar *line = new uchar[_n];
  memcpy(line, &matrix[i * _m], _n);
  return line;
}


BitMatrix *BitMatrix::transpose()
{
  uchar *mt = new uchar[_nb * _m];
  __sse_trans(matrix, mt, _nb, _mb);
  return new BitMatrix(mt, _mb, _n, !_le);
}


void BitMatrix::clear()
{
  memset(matrix, 0, _nb * _m);
}


void BitMatrix::dump(string file)
{
  if ( !matrix )
    throw runtime_error("Matrix is null");

  ofstream fout(file, ios::out | ios::binary);
  fout.write((char *) matrix, _nb * _m);
  fout.close();
}


void BitMatrix::print_bytes()
{
  cout << "\n\n";
  for ( size_t i = 0; i < _nb; i++ )
  {
    for ( size_t j = 0; j < _m; j++ )
      cout << "0x" << setfill('0') << setw(2) << hex << static_cast<int>(matrix[i * _m + j]) << " ";
    cout << "\n";
  }
  cout << endl;
}


void BitMatrix::print_bits()
{
  cout << "\n\n";
  string v;
  uchar tmp;
  for ( size_t i = 0; i < _nb; i++ )
  {
    for ( size_t j = 0; j < _m; j++ )
    {
      tmp = matrix[i * _m + j];
      for ( int p = 7; p >= 0; p-- )
      {
        v = ((tmp >> p) & 0x1) ? "1" : "0";
        cout << v;
      }
      cout << " ";
    }
    cout << "\n";
  }
  cout << endl;
}


// from https://mischasan.wordpress.com/2011/10/03/the-full-sse2-bit-matrix-transpose-routine/
void __sse_trans(uint8_t const *inp, uint8_t *out, int nrows, int ncols)
{
#   define INP(x, y) inp[(x)*ncols/8 + (y)/8]
#   define OUT(x, y) out[(y)*nrows/8 + (x)/8]
  ssize_t rr, cc, i, h;
  union
  {
    __m128i x;
    uint8_t b[16];
  } tmp;
  assert(nrows % 8 == 0 && ncols % 8 == 0);

#pragma omp parallel for private(rr, cc, i, tmp)
  // Do the main body in 16x8 blocks:
  for ( rr = 0; rr <= nrows - 16; rr += 16 )
  {
    for ( cc = 0; cc < ncols; cc += 8 )
    {
      for ( i = 0; i < 16; ++i )
        tmp.b[i] = INP(rr + i, cc);
      for ( i = 8; --i >= 0; tmp.x = _mm_slli_epi64(tmp.x, 1))
        *(uint16_t *) &OUT(rr, cc + i) = _mm_movemask_epi8(tmp.x);
    }
  }

  if ( nrows % 16 == 0 )
    return;
  rr = nrows - nrows % 16;

  // The remainder is a block of 8x(16n+8) bits (n may be 0).
  //  Do a PAIR of 8x8 blocks in each step:
  for ( cc = 0; cc <= ncols - 16; cc += 16 )
  {
    for ( i = 0; i < 8; ++i )
    {
      tmp.b[i] = h = *(uint16_t const *) &INP(rr + i, cc);
      tmp.b[i + 8] = h >> 8;
    }
    for ( i = 8; --i >= 0; tmp.x = _mm_slli_epi64(tmp.x, 1))
    {
      OUT(rr, cc + i) = h = _mm_movemask_epi8(tmp.x);
      OUT(rr, cc + i + 8) = h >> 8;
    }
  }
  if ( cc == ncols )
    return;

  //  Do the remaining 8x8 block:
  for ( i = 0; i < 8; ++i )
    tmp.b[i] = INP(rr + i, cc);
  for ( i = 8; --i >= 0; tmp.x = _mm_slli_epi64(tmp.x, 1))
    OUT(rr, cc + i) = _mm_movemask_epi8(tmp.x);
}
#endif
}; // end of namespace km