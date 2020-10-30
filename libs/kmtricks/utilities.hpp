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
#include <iostream>
#include <sstream>

#define NBYTE(bits) (((bits) >> 3) + ((bits) % 8 != 0))
#define NMOD8(byte) ((byte)+(8-((byte)%8)))

#define BITMASK(b) (1 << ((b) % CHAR_BIT))
#define BITSLOT(b) ((b) / CHAR_BIT)
#define BITSET(a, b) ((a)[BITSLOT(b)] |= BITMASK(b))

template <typename X>
void split(const std::string &s, char delim, X res)
{
  std::istringstream iss(s);
  std::string item;                                                                
  while (std::getline(iss, item, delim))
    *res++ = item;
}

std::vector<std::string> split(const std::string &s, char delim)
{
  std::vector<std::string> elems;
  split(s, delim, std::back_inserter(elems));
  return elems;
}

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
