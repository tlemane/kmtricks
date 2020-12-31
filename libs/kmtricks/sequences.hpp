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
 * 
*****************************************************************************/

#pragma once
#include <string>
#include <iostream>
#include <limits>
#include <cstring>
#include <vector>
#include <iterator>
#include <algorithm>
#include <iomanip>
#include <string>
#include "code.hpp"

using namespace std;

//! \defgroup Sequence

//! \namespace km
//! The kmtricks library namespace
namespace km
{

#define DEFAULT_MINIMIZER 1000000000

//! \ingroup Sequence
//! \class Hasher
//! \brief interface to build hasher for kmer 
//! \see Kmer::hash
template<typename K>
class Hasher
{
public:
  Hasher() { };
  virtual ~Hasher() {};
  
  //! \brief interface to build hasher for kmer
  //! \param data 
  //! \param seed 64-bit seed
  //! \return 64-bit hash value
  //! \see Kmer::hash
  virtual uint64_t operator()(K data, uint64_t seed) = 0;
};

//! \ingroup Sequence
//! \class XorHasher
//! \brief xor hashing for kmer
//! \see Kmer::hash
//! \see Hasher
template<typename K>
class XorHasher : public Hasher<K>
{
public:
  XorHasher() { };
  ~XorHasher() { };

  uint64_t operator()(K data, uint64_t seed) override
  {
    uint64_t hash = seed;
    uint64_t key = static_cast<uint64_t>(data);
    hash ^= (hash << 7) ^ key * (hash >> 3) ^ (~((hash << 11) + (key ^ (hash >> 5))));
    hash = (~hash) + (hash << 21);
    hash = hash ^ (hash >> 24);
    hash = (hash + (hash << 3)) + (hash << 8);
    hash = hash ^ (hash >> 14);
    hash = (hash + (hash << 2)) + (hash << 4);
    hash = hash ^ (hash >> 28);
    hash = hash + (hash << 31);
    return hash;
  }
};

//! \ingroup Sequence
//! \class Validator
//! \brief interface to build validator for kmers and minimizers
//! \see DefaultMinimizerValidator
//! \see Minimizer::set_kmer
//! \see Minimizer::set_superk
template<typename K>
class Validator
{
public:
  Validator() = default;
  virtual ~Validator() = default;
  virtual bool operator()(K value, size_t size) = 0;
};

//! \ingroup Sequence
//! \class DefaultMinimizerValidator
//! \brief Default validator for minimizer
//! \see Validator
//! \see Minimizer::set_kmer
//! \see Minimizer::set_superk
template<typename K>
class DefaultMinimizerValidator : public Validator<K>
{
public:
  DefaultMinimizerValidator() = default;
  ~DefaultMinimizerValidator() = default;
  bool operator()(K value, size_t size) override
  {
    K mask1 = numeric_limits<K>::max() >> (((sizeof(K)*8)-(size*2))+4);
    K mask01 = (K)0x5555555555555555;
    K mask00 = mask01 & mask1;
    value = ~((value) | (value >> 2));
    value = ((value >> 1) & value) & mask00;
    return !(value != 0);
  }
};

 //! \ingroup Sequence
 //! \class Kmer
 //! \brief A kmer representation
 //! 
 //! Only in {A, C, G, T}*
template<typename K>
class Kmer
{
public:
  //! \brief Constructor, empty kmer
  //! \param canonical :
  //! \param encoding :
  explicit Kmer(bool canonical, Code<K> *encoding = nullptr);

  //! \brief Constructor, from string
  //! \param kmer : 
  //! \param canonical :
  //! \param encoding :
  explicit Kmer(const string& kmer, bool canonical, Code<K> *encoding = nullptr);

  //! \brief Constructor, from unsigned int
  //! \param kmer :
  //! \param kmer_size :
  //! \param canonical :
  //! \param encoding :
  explicit Kmer(K kmer, size_t kmer_size, bool canonical, Code<K> *encoding = nullptr);

  //! \brief Destructor
  ~Kmer();

  //! \brief Copy Constructor
  Kmer(const Kmer<K> &k);

  //! \brief Copy Assignment Operator
  Kmer& operator=(const Kmer<K>&);

  //! \brief set kmer from string
  //! \param kmer : a kmer in string format, in {A,C,T,G}*
  void set_kmer(string kmer);
  
  //! \brief set kmer from unsigned int (2-bit encoded)
  //! \param kmer :
  //! \param kmer_size :
  void set_kmer(K kmer, size_t kmer_size);

  //! \brief get kmer int value
  //! \return kmer int value
  K value() const;
  
  //! \brief get kmer str value
  //! \return kmer as string
  string str_value();
  
  //! \brief get reverse complement ACCGA -> TCGGT
  //! \return reverse complement as int
  K rev_comp();

  //! \brief get reverse complement ACCGA -> TCGGT
  //! \return reverse complement as string
  string str_rev_comp();
  
  //! \brief check if the kmer is canonical
  //! \return true if canonical, else false
  bool is_canonical();

  //! \brief use canonical form
  void use_canonical();

  //! \brief get the kmer size
  //! \return kmer size
  uint size();

  //! \brief set default hasher (XorHasher)
  //! \see XorHasher
  //! \see Kmer::hash
  void set_default_hasher();

  //! \brief set a custom hasher
  //! \see Kmer::hash
  //! \see Hasher
  void set_hasher(Hasher<K> *hasher);

  //! \brief get kmer hash value, using internal hasher
  //! \param seed :
  //! \return hash value
  uint64_t hash(uint64_t seed = 0);

  //! \brief get kmer hash value, using external hasher
  //! \param hasher :
  //! \param seed :
  uint64_t hash(Hasher<K> *hasher, uint64_t seed = 0);

  //! \brief get kmer encoding
  //! \return Code
  //! \see Code
  Code<K> *get_encoding();

  //! \brief less than operator
  //! \param k
  //! \return true if this < that
  bool operator<(Kmer<K> const &k) const {return _bin_kmer < k.value();}
  //! \brief greater than operator
  //! \param k
  //! \return true if this > that
  bool operator>(Kmer<K> const &k) const {return _bin_kmer > k.value();}
  //! \brief not equal than operator
  //! \param k
  //! \return true != this != that
  bool operator!=(Kmer<K> const &k) const {return _bin_kmer != k.value();}
  //! \brief equal operator
  //! \param k
  //! \return true if this == that
  bool operator==(Kmer<K> const &k) const {return _bin_kmer == k.value();}

  
  //! \brief less than operator
  //! \param value
  //! \return true if this < that
  bool operator<(K const value) const {return _bin_kmer < value;}
  //! \brief greater than operator
  //! \param value
  //! \return true if this > that
  bool operator>(K const value) const {return _bin_kmer > value;}
  //! \brief not equal than operator
  //! \param value
  //! \return true != this != that
  bool operator!=(K const value) const {return _bin_kmer != value;}
  //! \brief equal operator
  //! \param value
  //! \return true if this == that
  bool operator==(K const value) const {return _bin_kmer == value;}


  //! \brief less than operator
  //! \param value
  //! \return true if this < that
  bool operator<(const string& value) const {return _bin_kmer < Kmer<K>(value, true).value();}
  //! \brief greater than operator
  //! \param value
  //! \return true if this > that
  bool operator>(const string& value) const {return _bin_kmer > Kmer<K>(value, true).value();}
  //! \brief not equal than operator
  //! \param value
  //! \return true != this != that
  bool operator!=(const string& value) const {return _bin_kmer != Kmer<K>(value, true).value();}
  //! \brief equal operator
  //! \param value
  //! \return true if this == that
  bool operator==(const string& value) const {return _bin_kmer == Kmer<K>(value, true).value();}

  ////! \brief less than operator
  ////! \param value
  ////! \return true if this < that
  //bool operator<(const char *value) const {return _bin_kmer < Kmer<K>(value, true).value();}
  ////! \brief greater than operator
  ////! \param value
  ////! \return true if this > that
  //bool operator>(const char *value) const {return _bin_kmer > Kmer<K>(value, true).value();}
  ////! \brief not equal than operator
  ////! \param value
  ////! \return true != this != that
  //bool operator!=(const char *value) const {return _bin_kmer != Kmer<K>(value, true).value();}
  ////! \brief equal operator
  ////! \param value
  ////! \return true if this == that
  //bool operator==(const char *value) const {return _bin_kmer == Kmer<K>(value, true).value();}

private:
  Code<K> *_code;
  bool    _custom_hash;

  bool    _has_bin;
  bool    _is_canonical;

  K       _bin_kmer;
  size_t  _size;
  K       _kmer_mask;

  Hasher<K> *_hasher;

  bool      _custom_enc;
};

//! \ingroup Sequence
//! \class Superk
//! \brief A superkmer representation
//!
//! Corresponds to a set of overlapping kmers that shared the same minimizer. \see Minimizer
template<typename K>
class Superk
{
public:
  //! \brief Constructor, empty superk
  //! \param kmer_size :
  //! \param encoding :
  explicit Superk(size_t kmer_size, Code<K> *encoding = nullptr);
  
  //! \brief Constructor, from string
  //! \param superkmer :
  //! \param kmer_size :
  //! \param encoding :
  explicit Superk(string superkmer, size_t kmer_size, Code<K> *encoding = nullptr);
  
  //! \brief Constructor, from buffer, 2-bit encoded array
  //! \param buffer :
  //! \param superk_size :
  //! \param kmer_size :
  //! \param gatb_format :
  //! \param encoding :
  explicit Superk(uchar* buffer, size_t superk_size, size_t kmer_size,
         bool gatb_format = false, Code<K> *encoding = nullptr);
  
  //! \brief Destructor
  ~Superk();

  //! \brief Copy Constructor
  Superk(const Superk<K> &s);

  //! \brief Copy Assignment Operator
  Superk<K>& operator=(const Superk<K> &s);

  string str_value() const;
  uchar  *value();
  Kmer<K> get_kmer(bool canonical);
  void    get_kmer(Kmer<K> *kmer);
  Kmer<K> get_kmer(int n, bool canonical);
  void    get_kmer(int n, Kmer<K> *kmer);

  Kmer<K> get_first();
  size_t  size();
  size_t  nb_kmers();

  Code<K> *get_encoding();

  void set_superk(string superkmer);
  void set_superk(uchar *buffer, size_t superk_size, size_t kmer_size, bool gatb_format = false);

  bool operator<(Superk<K> const &s) const {return str_value().compare(s.str_value()) < 0;}
  bool operator>(Superk<K> const &s) const {return str_value().compare(s.str_value()) > 0;}
  bool operator!=(Superk<K> const &s) const {return str_value().compare(s.str_value()) != 0;}
  bool operator==(Superk<K> const &s) const {return str_value().compare(s.str_value()) == 0;}

  //bool operator<(const char *value) const {return str_value().compare(value) < 0;}
  //bool operator>(const char *value) const {return str_value().compare(value) > 0;}
  //bool operator!=(const char *value) const {return str_value().compare(value) != 0;}
  //bool operator==(const char *value) const {return str_value().compare(value) == 0;}

  bool operator<(const string& value) const {return str_value().compare(value) < 0;}
  bool operator>(const string& value) const {return str_value().compare(value) > 0;}
  bool operator!=(const string& value) const {return str_value().compare(value) != 0;}
  bool operator==(const string& value) const {return str_value().compare(value) == 0;}

  void operator++(int) {_kmer_index++;};
  void operator--(int) {_kmer_index--;}

private:
  void _build_from_string(string superkmer);
  void _build_from_gatb_format(uchar *buffer);


private:
  Code<K> *_code;
  uchar   *_bin_superk;

  size_t  _ksize;
  size_t  _superksize;

  size_t  _kmer_index;

  K       _kmer_mask;
  bool    _gatb;

  bool    _custom_enc;
};

//! \ingroup Sequence
//! \class Minimizer
//! \brief A minimizer representation
//! 
//! For a sequence S, corresponds to the lexicographically smaller sub-sequence of size m
template<typename K>
class Minimizer
{
public:
  //! \brief Constructor, empty minimizer
  //! \param msize :
  //! \param default_minim :
  //! \param validator :
  Minimizer(size_t msize, K default_minim = DEFAULT_MINIMIZER, Validator<K> *validator = nullptr);

  //! \brief Constructor, from kmer with default validator
  //! \see Kmer
  //! \param kmer :
  //! \param msize :
  //! \param check_validity :
  //! \param default_minim :
  Minimizer(Kmer<K> *kmer, size_t msize, bool check_validity, K default_minim = DEFAULT_MINIMIZER);
  
  //! \brief Constructor, from kmer with custom validator
  //! \see Kmer
  //! \param kmer :
  //! \param msize :
  //! \param validator :
  //! \param default_minim :
  Minimizer(Kmer<K> *kmer, size_t msize, Validator<K> *validator = nullptr, K default_minim = DEFAULT_MINIMIZER);

  //! \brief Constructor, from superk with default validator
  //! \see Superk
  //! \param superk :
  //! \param msize :
  //! \param check_validity :
  //! \param default_minim :
  Minimizer(Superk<K> *superk, size_t msize, bool check_validity, K default_minim = DEFAULT_MINIMIZER);

  //! \brief Constructor, from superk with custom validator
  //! \see Superk
  //! \param superk
  //! \param msize :
  //! \param validator :
  //! \param default_minim :
  Minimizer(Superk<K> *superk, size_t msize, Validator<K> *validator = nullptr, K default_minim = DEFAULT_MINIMIZER);

  //! \brief Destructor
  ~Minimizer();

  //! \brief Copy Constructor
  Minimizer(const Minimizer<K> &m);

  //! \brief Copy Assignement Operator
  Minimizer<K>& operator=(const Minimizer<K> &m);

  K       value() const;
  string  str_value();

  void set_default();
  void set_default(K minimizer);
  void set_default(string minimizer);
  void set_kmer(Kmer<K> *kmer, size_t msize, bool check_validity);
  void set_superk(Superk<K> *superk, size_t msize, bool check_validity);

  bool operator<(Minimizer<K> const &m) const {return _minimizer < m.value();}
  bool operator>(Minimizer<K> const &m) const {return _minimizer > m.value();}
  bool operator!=(Minimizer<K> const &m) const {return _minimizer != m.value();}
  bool operator==(Minimizer<K> const &m) const {return _minimizer == m.value();}

  bool operator<(K const value) const {return _minimizer < value;}
  bool operator>(K const value) const {return _minimizer > value;}
  bool operator!=(K const value) const {return _minimizer != value;}
  bool operator==(K const value) const {return _minimizer == value;}

  bool operator<(const string& value) const
  {return _minimizer < _code->encode(value, value.size());}
  bool operator>(const string& value) const
  {return _minimizer > _code->encode(value, value.size());}
  bool operator!=(const string& value) const
  {return _minimizer != _code->encode(value, value.size());}
  bool operator==(const string& value) const
  {return _minimizer == _code->encode(value, value.size());}

  //bool operator<(const char *value) const {return _minimizer < _code->encode(value, strlen(value));}
  //bool operator>(const char *value) const {return _minimizer > _code->encode(value, strlen(value));}
  //bool operator!=(const char *value) const {return _minimizer != _code->encode(value, strlen(value));}
  //bool operator==(const char *value) const {return _minimizer == _code->encode(value, strlen(value));}

private:
  void _minimizer_from_kmer();
  void _minimizer_from_superk();
  K    _rev_comp(const K& minimizer);

private:
  Kmer<K>     *_kmer;
  Superk<K>   *_superk;
  Code<K>     *_code;
  size_t      _msize;
  K           _minimizer;
  bool        _check;
  K           _default;

  Validator<K> *_is_valid_minimizer;
  bool        _out_valid;
};

#ifndef _KM_LIB_INCLUDE_
///////////////
// KMER IMPL //
///////////////

template<typename K>
Kmer<K>::Kmer(bool canonical, Code<K> *encoding)
  :
  _code(encoding),
  _custom_hash(false),
  _has_bin(false),
  _is_canonical(canonical),
  _bin_kmer(0),
  _size(0),
  _kmer_mask(0),
  _hasher(nullptr)
{
  if (!_code)
  {
    _code = new Code<K>();
    _custom_enc = false;
  }
  else
  {
    _code = encoding;
    _custom_enc = true;
  }
  set_default_hasher();
}


template<typename K>
Kmer<K>::Kmer(const string& kmer, bool canonical, Code<K> *encoding)
  :
  _code(encoding),
  _custom_hash(false),
  _has_bin(false),
  _is_canonical(canonical),
  _bin_kmer(0),
  _size(0),
  _kmer_mask(0),
  _hasher(nullptr)
{
  if (!_code)
  {
    _code = new Code<K>();
    _custom_enc = false;
  }
  else
  {
    _code = encoding;
    _custom_enc = true;
  }
  set_default_hasher();
  set_kmer(kmer);
}


template<typename K>
Kmer<K>::Kmer(K kmer, size_t kmer_size, bool canonical, Code<K> *encoding)
  :
  _code(encoding),
  _custom_hash(false),
  _has_bin(false),
  _is_canonical(canonical),
  _bin_kmer(0),
  _size(kmer_size),
  _kmer_mask(0),
  _hasher(nullptr)
{
  if (!_code)
  {
    _code = new Code<K>();
    _custom_enc = false;
  }
  else
  {
    _code = encoding;
    _custom_enc = true;
  }
  set_default_hasher();
  set_kmer(kmer, kmer_size);
}


template<typename K>
Kmer<K>::Kmer(const Kmer<K> &k)
  : _code(k._code),
    _custom_hash(k._custom_hash),
    _has_bin(k._has_bin),
    _is_canonical(k._is_canonical),
    _bin_kmer(k._bin_kmer),
    _size(k._size),
    _kmer_mask(k._kmer_mask),
    _hasher(k._hasher),
    _custom_enc(k._custom_enc)
{
if (!_custom_enc)
    _code = new Code<K>();
  if (!_custom_hash)
    _hasher = new XorHasher<K>();
}

template<typename K>
Kmer<K>& Kmer<K>::operator=(const Kmer<K> &k)
{
  if (this != &k)
  {
    _code         = k._code;
    _custom_hash  = k._custom_hash;
    _custom_enc   = k._custom_enc;
    _has_bin      = k._has_bin;
    _is_canonical = k._is_canonical;
    _bin_kmer     = k._bin_kmer;
    _size         = k._size;
    _kmer_mask    = k._kmer_mask;
    _hasher       = k._hasher;

    if (!_custom_enc)
      _code = new Code<K>();
    if (!_custom_hash)
      _hasher = new XorHasher<K>();
  }
  return *this;
}

template<typename K>
Kmer<K>::~Kmer()
{
  if (!_custom_enc)  delete _code;
  if (!_custom_hash) delete _hasher;
}


template<typename K>
uint Kmer<K>::size()
{
  return _size;
}


template<typename K>
K Kmer<K>::value() const
{
  return _bin_kmer;
}


template<typename K>
string Kmer<K>::str_value()
{
  return _code->decode(_bin_kmer, _size);
}


template<typename K>
bool Kmer<K>::is_canonical()
{
  return _is_canonical;
}


template<typename K>
Code<K>* Kmer<K>::get_encoding()
{
  return _code;
}


template<typename K>
void Kmer<K>::set_default_hasher()
{
  if ( !_custom_hash )
    if (_hasher)
      delete _hasher;
  _hasher = new XorHasher<K>();
  _custom_hash = false;
}


template<typename K>
void Kmer<K>::set_hasher(Hasher<K> *hasher)
{
  if (!_custom_hash)
    if ( _hasher )
      delete _hasher;
  _hasher = hasher;
  _custom_hash = true;
}


template<typename K>
void Kmer<K>::set_kmer(string kmer)
{
  _kmer_mask = numeric_limits<K>::max() >> ((sizeof(K) * 8) - (_size * 2));
  _size = kmer.size();
  _bin_kmer = _code->encode(kmer, _size);
  _has_bin = true;

  if ( _is_canonical )
  {
    K rev = rev_comp();
    if ( rev < _bin_kmer )
      _bin_kmer = rev;
  }
}


template<typename K>
void Kmer<K>::set_kmer(K kmer, size_t kmer_size)
{
  _kmer_mask = numeric_limits<K>::max() >> ((sizeof(K) * 8) - (_size * 2));
  _size = kmer_size;
  _bin_kmer = kmer & _kmer_mask;
  _has_bin = true;

  if ( _is_canonical )
  {
    K rev = rev_comp();
    if ( rev < _bin_kmer )
      _bin_kmer = rev;
  }
}


template<typename K>
void Kmer<K>::use_canonical()
{
  _is_canonical = true;
  if ( _has_bin )
  {
    K rev = rev_comp();
    if ( rev < _bin_kmer )
      _bin_kmer = rev;
  }
}


template<typename K>
uint64_t Kmer<K>::hash(uint64_t seed)
{
  return (*_hasher)(_bin_kmer, seed);
}


template<typename K>
uint64_t Kmer<K>::hash(Hasher<K> *hasher, uint64_t seed)
{
  return (*hasher)(_bin_kmer, seed);
}


template<typename K>
K Kmer<K>::rev_comp()
{
  K res = 0;
  K seq = _bin_kmer;
  for ( int c = _size - 1; c >= 0; c-- )
  {
    res <<= 2;
    res |= _code->_NToB[_code->_revC[seq & 3]];
    seq >>= 2;
  }
  return res & _kmer_mask;
}


template<typename K>
string Kmer<K>::str_rev_comp()
{
  return _code->decode(rev_comp(), _size);
}


/////////////////
// SUPERK IMPL //
/////////////////

template<typename K>
Superk<K>::Superk(size_t kmer_size, Code<K> *encoding)
  :
  _code(encoding),
  _ksize(kmer_size),
  _superksize(0),
  _kmer_index(0)
{
  _bin_superk = static_cast<uchar*>(calloc(1, sizeof(uchar)));
  if (!_code)
  {
    _code = new Code<K>();
    _custom_enc = false;
  }
  else
  {
    _code = encoding;
    _custom_enc = true;
  }
  _kmer_mask = numeric_limits<K>::max() >> ((sizeof(K) * 8) - (kmer_size * 2));
}


template<typename K>
Superk<K>::Superk(string superkmer, size_t kmer_size, Code<K> *encoding)
  :
  _code(encoding),
  _ksize(kmer_size),
  _kmer_index(0)
{
  _superksize = superkmer.size();
  _bin_superk = static_cast<uchar*>(calloc((_superksize/4)+1, sizeof(uchar)));
  _kmer_mask = numeric_limits<K>::max() >> ((sizeof(K) * 8) - (kmer_size * 2));
  if (!_code)
  {
    _code = new Code<K>();
    _custom_enc = false;
  }
  else
  {
    _code = encoding;
    _custom_enc = true;
  }
  _build_from_string(superkmer);
}


template<typename K>
Superk<K>::Superk(uchar* buffer, size_t superk_size, size_t kmer_size, bool gatb_format, Code<K> *encoding)
  :
  _code(encoding),
  _bin_superk(nullptr),
  _superksize(superk_size),
  _kmer_index(0),
  _ksize(kmer_size),
  _gatb(gatb_format)
{
  int nb_bytes = (_superksize/4)+1;
  _bin_superk = static_cast<uchar*>(calloc(nb_bytes, sizeof(uchar)));
  _kmer_mask = numeric_limits<K>::max() >> ((sizeof(K) * 8) - (kmer_size * 2));
  if (!_code)
  {
    _code = new Code<K>;
    _custom_enc = false;
  }
  else
  {
    _code = encoding;
    _custom_enc = true;
  }
  if (!_gatb)
    memcpy(_bin_superk, buffer, nb_bytes);
  else
    _build_from_gatb_format(buffer);
}


template<typename K>
Superk<K>::~Superk()
{
  free (_bin_superk);
  if (!_custom_enc) delete _code;
}


template<typename K>
Superk<K>::Superk(const Superk<K> &s)
: _code(nullptr),
  _bin_superk(nullptr),
  _superksize(s._superksize),
  _kmer_index(s._kmer_index),
  _ksize(s._ksize),
  _gatb(s._gatb),
  _custom_enc(s._custom_enc),
  _kmer_mask(s._kmer_mask)
{
  if (!_custom_enc)
    _code = new Code<K>();

  int nb_bytes = (_superksize/4)+1;
  _bin_superk = static_cast<uchar*>(calloc(nb_bytes, sizeof(uchar)));
  memcpy(_bin_superk, s._bin_superk, nb_bytes);
}


template<typename K>
Superk<K>& Superk<K>::operator=(const Superk<K> &s)
{
  _code = nullptr;
  _bin_superk = nullptr;
  _superksize = s._superksize;
  _kmer_index = s._kmer_index;
  _ksize = s._ksize;
  _gatb = s._gatb;
  _custom_enc = s._custom_enc;
  _kmer_mask = s._kmer_mask;

  if (!_custom_enc)
    _code = new Code<K>();
  int nb_bytes = (_superksize/4)+1;
  _bin_superk = static_cast<uchar*>(calloc(nb_bytes, sizeof(uchar)));
  memcpy(_bin_superk, s._bin_superk, nb_bytes);
  return *this;
}


template<typename K>
void Superk<K>::_build_from_string(string superkmer)
{
  for (size_t i=0; i<_superksize; i++)
  {
    int pos = i/4;
    _bin_superk[pos] <<= 2;
    _bin_superk[pos] |= static_cast<uchar>(_code->encode(superkmer[i]));
  }
  int shift = (_superksize%4);
  _bin_superk[(_superksize-1)/4] <<=  shift ? 8-(_superksize%4)*2 : 0;
}


template<typename K>
void Superk<K>::_build_from_gatb_format(uchar *buffer)
{
  uchar *ptr = buffer;
  uint8_t nbK = _superksize-_ksize+1;
  uchar newbyte = 0;

  int rem_size = _ksize;

  K first_kmer = 0;
  K Tnewbyte;

  int nbr = 0;
  int index = _ksize % 4 ? (_ksize/4) : (_ksize/4)-1;

  while (rem_size >= 4)
  {
    newbyte = *ptr; ptr++;
    _bin_superk[index] = newbyte; index--;
    Tnewbyte = newbyte;
    first_kmer = first_kmer | (Tnewbyte << (8*nbr));
    rem_size -= 4; nbr++;
  }

  int uid = 4;

  if (rem_size>0)
  {
    newbyte = *ptr ; ptr++;
    _bin_superk[index] = newbyte; index--;
    Tnewbyte = newbyte;
    first_kmer = first_kmer | (Tnewbyte << (8*nbr));
    uid = rem_size;
  }

  first_kmer = first_kmer & _kmer_mask;
  uint8_t rem = nbK;
  uchar newnt = 0;

  size_t curr_offset = _ksize/4;
  size_t shift_size = _ksize % 4 ? (6-(2*((_ksize%4)-1))) : 6 ;

  _bin_superk[curr_offset] >>= shift_size;
  int nbnt = _ksize%4;
  for (int i=0; i<nbK; i++, rem--)
  {
    if (rem < 2) break;

    if (uid >= 4)
    {
      newbyte = *ptr; ptr++;
      Tnewbyte = newbyte;
      uid = 0;
    }
    newnt = (Tnewbyte >> (2*uid)) & 3; uid++;
    _bin_superk[curr_offset] <<= 2;
    _bin_superk[curr_offset] |= newnt;
    nbnt++;
    if (nbnt == 4)
    {
      curr_offset++;
      nbnt = 0;
    }
  }
  while (nbnt<4)
  {
    _bin_superk[curr_offset] <<= 2;
    nbnt++;
  }
}


template<typename K>
void Superk<K>::set_superk(string superkmer)
{
  size_t new_size = superkmer.size();
  if (new_size > _superksize)
    _bin_superk = static_cast<uchar*>(realloc(_bin_superk, (new_size/4)+1));
  _superksize = new_size;
  _build_from_string(superkmer);
}


template<typename K>
void Superk<K>::set_superk(uchar *buffer, size_t superk_size, size_t kmer_size, bool gatb_format)
{
  _ksize = kmer_size;
  int nb_bytes = (superk_size/4)+1;
  if (superk_size > _superksize)
    _bin_superk = static_cast<uchar*>(realloc(_bin_superk, nb_bytes));
  _superksize = superk_size;
  if (!gatb_format)
    memcpy(_bin_superk, buffer, nb_bytes);
  else
    _build_from_gatb_format(buffer);
}


template<typename K>
uchar* Superk<K>::value()
{
  return _bin_superk;
}


template<typename K>
string Superk<K>::str_value() const
{
  string ret;
  ret = "";
  size_t nb_bytes = (_superksize/4)+1;
  for (size_t i=0; i<nb_bytes; i++)
  {
    ret += _code->decode(_bin_superk[i]);
  }
  return ret.substr(0, _superksize);
}


template<typename K>
size_t Superk<K>::size()
{
  return _superksize;
}


template<typename K>
size_t Superk<K>::nb_kmers()
{
  return _superksize-_ksize+1;
}


template<typename K>
Code<K> *Superk<K>::get_encoding()
{
  return _code;
}


template<typename K>
Kmer<K> Superk<K>::get_first()
{
  size_t tmp_size = 0;
  size_t end      = (_ksize-1)/4;

  K value = 0;

  for (size_t i=0; i<=end; i++)
  {
    for (int j=6; j>=0; j-=2)
    {
      value <<= 2;
      value |= ((_bin_superk[i] >> j) & 3);
      tmp_size++;
      if (tmp_size == _ksize)
      {
        tmp_size = 0;
        break;
      }
    }
    if (!tmp_size) break;
  }

  return Kmer<K>(value, _ksize, true, _code);
}


template<typename K>
Kmer<K> Superk<K>::get_kmer(int n, bool canonical)
{
  size_t tmp_size = 0;
  K value = 0;
  int start = n/4;
  int end = start+(_ksize/4)+1;

  bool first = true;
  for ( int i = start; i <= end; i++ )
  {
    int s = 6;
    if (first)
    {
      s = 6-(2*(n%4));
      first = false;
    }

    for ( int j = s; j >= 0; j -= 2 )
    {
      value <<= 2;
      value |= ((_bin_superk[i] >> j) & 3);
      tmp_size++;

      if (tmp_size == _ksize)
      {
        tmp_size = 0;
        break;
      }
    }
    if (tmp_size == 0)
      break;
  }
  return Kmer<K>(value, _ksize, canonical, _code);
}


template<typename K>
void Superk<K>::get_kmer(int n, Kmer<K> *kmer)
{
  size_t tmp_size = 0;
  K value = 0;
  int start = n/4;
  int end = start+(_ksize/4)+1;

  bool first = true;
  for ( int i = start; i <= end; i++ )
  {
    int s = 6;
    if (first)
    {
      s = 6-(2*(n%4));
      first = false;
    }

    for ( int j = s; j >= 0; j -= 2 )
    {
      value <<= 2;
      value |= ((_bin_superk[i] >> j) & 3);
      tmp_size++;

      if (tmp_size == _ksize)
      {
        tmp_size = 0;
        break;
      }
    }
    if (tmp_size == 0)
    {
      kmer->set_kmer(value, _ksize);
    }
  }
}


template<typename K>
Kmer<K> Superk<K>::get_kmer(bool canonical)
{
  return get_kmer(_kmer_index, canonical);
}


template<typename K>
void Superk<K>::get_kmer(Kmer<K> *kmer)
{
  get_kmer(_kmer_index, kmer);
}
////////////////////
// MINIMIZER IMPL //
////////////////////

template<typename K>
Minimizer<K>::~Minimizer()
{
  if (_is_valid_minimizer && !_out_valid)
    delete _is_valid_minimizer;
}


template<typename K>
Minimizer<K>::Minimizer(size_t msize, K default_minim, Validator<K> *validator)
: _kmer(nullptr),
  _superk(nullptr),
  _code(nullptr),
  _msize(msize),
  _default(default_minim),
  _is_valid_minimizer(validator),
  _out_valid(false)
{
  _minimizer = numeric_limits<K>::max();
  if (!_is_valid_minimizer) _is_valid_minimizer = new DefaultMinimizerValidator<K>();
  else _out_valid  = true;
}


template<typename K>
Minimizer<K>::Minimizer(Kmer<K> *kmer, size_t msize, bool check_validity, K default_minim)
: _kmer(kmer),
  _superk(nullptr),
  _msize(msize),
  _check(check_validity),
  _default(default_minim),
  _is_valid_minimizer(nullptr),
  _out_valid(false)
{
  _code = _kmer->get_encoding();
  _minimizer = numeric_limits<K>::max();
  if (_check) _is_valid_minimizer = new DefaultMinimizerValidator<K>();
  _minimizer_from_kmer();
}


template<typename K>
Minimizer<K>::Minimizer(Kmer<K> *kmer, size_t msize, Validator<K> *validator, K default_minim)
: _kmer(kmer),
  _superk(nullptr),
  _msize(msize),
  _check(true),
  _default(default_minim),
  _is_valid_minimizer(validator),
  _out_valid(false)
{
  _code = _kmer->get_encoding();
  _minimizer = numeric_limits<K>::max();
  if (!_is_valid_minimizer) _is_valid_minimizer = new DefaultMinimizerValidator<K>();
  else _out_valid = true;
  _minimizer_from_kmer();
}


template<typename K>
Minimizer<K>::Minimizer(Superk<K> *superk, size_t msize, bool check_validity, K default_minim)
: _kmer(nullptr),
  _superk(superk),
  _msize(msize),
  _check(check_validity),
  _default(default_minim),
  _is_valid_minimizer(nullptr),
  _out_valid(false)
{
  _code = _superk->get_encoding();
  _minimizer = numeric_limits<K>::max();
  if (_check) _is_valid_minimizer = new DefaultMinimizerValidator<K>();
  _minimizer_from_superk();
}


template<typename K>
Minimizer<K>::Minimizer(Superk<K> *superk, size_t msize, Validator<K> *validator, K default_minim)
: _kmer(nullptr),
  _superk(superk),
  _msize(msize),
  _default(default_minim),
  _is_valid_minimizer(validator),
  _check(true),
  _out_valid(false)
{
  _code = _superk->get_encoding();
  _minimizer = numeric_limits<K>::max();
  if (!_is_valid_minimizer) _is_valid_minimizer = new DefaultMinimizerValidator<K>();
  else _out_valid = true;
  _minimizer_from_superk();
}


template<typename K>
Minimizer<K>::Minimizer(const Minimizer<K> &m)
: _kmer(m._kmer),
  _superk(m._superk),
  _msize(m._msize),
  _default(m._default),
  _is_valid_minimizer(m._is_valid_minimizer),
  _check(m._check),
  _out_valid(m._out_valid)
{
  if (_is_valid_minimizer && !_out_valid)
    _is_valid_minimizer = new DefaultMinimizerValidator<K>();
}


template<typename K>
Minimizer<K>& Minimizer<K>::operator=(const Minimizer<K> &m)
{
  _kmer = m._kmer;
  _superk = m._superk;
  _msize = m._msize;
  _default = m._default;
  _is_valid_minimizer = m._is_valid_minimizer;
  _check = m._check;
  _out_valid = m._out_valid;

  if (_is_valid_minimizer && !_out_valid)
    _is_valid_minimizer = new DefaultMinimizerValidator<K>();
  return *this;
}


template<typename K>
void Minimizer<K>::set_default()
{
  _default = DEFAULT_MINIMIZER;
  _minimizer = numeric_limits<K>::max();
  if (_kmer)
    _minimizer_from_kmer();
  else if (_superk)
    _minimizer_from_superk();
}


template<typename K>
void Minimizer<K>::set_default(K minimizer)
{
  _default = minimizer;
  _minimizer = numeric_limits<K>::max();
  if (_kmer)
    _minimizer_from_kmer();
  else if (_superk)
    _minimizer_from_superk();
}


template<typename K>
void Minimizer<K>::set_default(string minimizer)
{
  size_t msize = minimizer.size();
  if (msize != _msize)
    throw runtime_error("Invalid minimizer size, str size = " + to_string(msize) + ".");
  _default = _code->encode(minimizer, msize);
  _minimizer = numeric_limits<K>::max();

  if (_kmer)
    _minimizer_from_kmer();
  else if (_superk)
  {
    _minimizer_from_superk();
  }
}

template<typename K>
K Minimizer<K>::_rev_comp(const K& minimizer)
{
  K res = 0;
  K seq = minimizer;
  K _mmer_mask = numeric_limits<K>::max() >> ((sizeof(K)*8)-(_msize*2));
  for ( int c = _msize - 1; c >= 0; c-- )
  {
    res <<= 2;
    res |= _code->_NToB[_code->_revC[seq & 3]];
    seq >>= 2;
  }
  return res & _mmer_mask;
}

template<typename K>
void Minimizer<K>::_minimizer_from_kmer()
{
  if (!_kmer->is_canonical())
    cerr << "Warning: minimzer on non-canonical k-mer" << endl;

  K _mmer_mask = numeric_limits<K>::max() >> ((sizeof(K)*8)-(_msize*2));
  K _bink = _kmer->value();
  int nb_mmers = _kmer->size() - _msize + 1;
  for (int i=nb_mmers-1; i>=0; i--)
  {
    K tmp = (_bink >> (i*2)) & _mmer_mask;
    K rev = _rev_comp(tmp);
    if (rev < tmp) tmp = rev;
    if (_check)
    {
      if ((*_is_valid_minimizer)(tmp, _msize) && tmp < _minimizer)
        _minimizer = tmp;
    }
    else
      if (tmp < _minimizer) _minimizer = tmp;
  }
  if (_minimizer == numeric_limits<K>::max())
    _minimizer = _default;
}


template<typename K>
void Minimizer<K>::_minimizer_from_superk()
{
  Kmer<K> kmer = _superk->get_first();
  if (!kmer.is_canonical())
    cerr << "Warning: minimzer on non-canonical k-mer" << endl;

  K _mmer_mask = numeric_limits<K>::max() >> ((sizeof(K)*8)-(_msize*2));
  K _bink = kmer.value();
  int nb_mmers = kmer.size() - _msize + 1;

  for (int i=nb_mmers-1; i>=0; i--)
  {
    K tmp = (_bink >> (i*2)) & _mmer_mask;
    if (_check)
    {
      if ((*_is_valid_minimizer)(tmp, _msize) && tmp < _minimizer)
        _minimizer = tmp;
    }
    else
    if (tmp < _minimizer) _minimizer = tmp;
  }

  if (_minimizer == numeric_limits<K>::max())
    _minimizer = _default;
}


template<typename K>
void Minimizer<K>::set_kmer(Kmer<K> *kmer, size_t msize, bool check_validity)
{
  _kmer = kmer;
  _code = kmer->get_encoding();
  _msize = msize;
  _check = check_validity;
  _minimizer = numeric_limits<K>::max();
  _minimizer_from_kmer();
}

template<typename K>
void Minimizer<K>::set_superk(Superk<K> *superk, size_t msize, bool check_validity)
{
  _superk = superk;
  _code = superk->get_encoding();
  _msize = msize;
  _check = check_validity;
  _minimizer = numeric_limits<K>::max();
  _minimizer_from_superk();
}

template<typename K>
string Minimizer<K>::str_value()
{
  return _code->decode(_minimizer, _msize);
}


template<typename K>
K Minimizer<K>::value() const
{
  return _minimizer;
}
#endif
} // end of namespace km