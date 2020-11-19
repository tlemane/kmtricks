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
#include <algorithm>
#include <fstream>
#include <tuple>
#include <iterator>
#include <stdexcept>
#include <exception>

/*! \defgroup Utilities
*/

//!\ingroup Utilities
#define NBYTE(bits) (((bits) >> 3) + ((bits) % 8 != 0))
//!\ingroup Utilities
#define NMOD8(byte) ((byte)+(8-((byte)%8)))

//!\ingroup Utilities
#define BITMASK(b) (1 << ((b) % CHAR_BIT))
//!\ingroup Utilities
#define BITSLOT(b) ((b) / CHAR_BIT)
//!\ingroup Utilities
#define BITSET(a, b) ((a)[BITSLOT(b)] |= BITMASK(b))

template <typename X>
void split(const std::string &s, char delim, X res)
{
  std::istringstream iss(s);
  std::string item;                                                                
  while (std::getline(iss, item, delim))
    *res++ = item;
}

template<typename Callback>
std::vector<std::string> split_strip(const std::string&s, char delim, Callback func)
{
  std::vector<std::string> elems;
  std::istringstream iss(s);
  std::string item;
  auto res = std::back_inserter(elems);
  while (std::getline(iss, item, delim))
    *res++ = func(item);
  return elems;
}

/*! \ingroup Utilities 
*   \brief split string into std::vector<std::string>
*/
std::vector<std::string> split(const std::string &s, char delim)
{
  std::vector<std::string> elems;
  split(s, delim, std::back_inserter(elems));
  return elems;
}

using fof_t = std::vector<std::tuple<std::string, std::vector<std::string>, int>>;

fof_t parse_km_fof(std::string fof_path)
{
  fof_t fof_vec;
  std::ifstream in_fof(fof_path);
  if (!in_fof.good()) throw std::runtime_error("Unable to open fof at " + fof_path + ".");
  std::string line;
  
  while (getline(in_fof, line) && line.size() > 0)
  {
    auto rspace = [] (std::string s) -> std::string {
      s.erase(remove(s.begin(), s.end(), ' '), s.end()); 
      return s;
    };
    
    if (rspace(line).size() == 0) continue;
    std::vector<std::string> temp = split_strip(line, ':', rspace);
    if (temp.size() != 2) throw std::runtime_error("Fof bad format1.");
    std::string id = temp[0];
    id.erase(remove(id.begin(), id.end(), ' '), id.end());

    std::vector<std::string> temp2 = split_strip(temp[1], '!', rspace);
    if (temp2.size() != 1 && temp2.size() != 2) throw std::runtime_error("Fof bad format2.");
    int count = -1;
    if (temp2.size() == 2) count = stoi(temp2[1]);
    
    std::vector<std::string> files = split_strip(temp2[0], ';', rspace);
    fof_vec.push_back(make_tuple(id, files, count));
  }
  return fof_vec;
}

std::string all_files(fof_t& fof)
{
  std::vector<std::string> allv;
  for (auto& elem: fof)
    for (auto& f: std::get<1>(elem))
      allv.push_back(f);
  std::ostringstream alls;
  const char* delim = ",";
  std::copy(allv.begin(), allv.end(), std::ostream_iterator<std::string>(alls, delim));
  return alls.str();
}

#ifdef __SIZEOF_INT128__
//! \ingroup Utilities
//! \brief __uint128_t support for std::ostream
std::ostream&
operator<<( std::ostream& dest, __uint128_t value )
{
    const char digit[11] = "0123456789";
    std::ostream::sentry s( dest );
    if ( s ) {
        char buffer[ 128 ];
        char* d = std::end( buffer );
        do
        {
            --d;
            *d = digit[ value % 10 ];
            value/= 10;
        } while ( value != 0 );
        int len = std::end( buffer ) - d;
        if ( dest.rdbuf()->sputn( d, len ) != len ) {
            dest.setstate( std::ios_base::badbit );
        }
    }
    return dest;
}
#endif

//! \cond HIDDEN_SYMBOLS
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

//! \endcond

//! \ingroup Utilities 
//! \brief select kmer int type from k
//!
//! typedef selectK<31>::type   KType;
template<unsigned long long klength>
struct selectK : select_<requiredK<klength>::value> {};

//! \ingroup Utilities 
//! \brief select count int type from max count
//!
//! typedef selectC<255>::type   CType;
template<unsigned long long klength>
struct selectC : select_<requiredC<klength>::value> {};
