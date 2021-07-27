/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
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

/** \file Integer.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Entry point class for large integer usage
 */
#ifndef _GATB_CORE_TOOLS_MATH_INTEGER_HPP_
#define _GATB_CORE_TOOLS_MATH_INTEGER_HPP_

/********************************************************************************/
#include <gatb/tools/math/LargeInt.hpp>
#include <gatb/system/api/Exception.hpp>

#include <boost/variant.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/pop_front.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/at.hpp>

/********************************************************************************/
namespace gatb  {  namespace core  { namespace tools {  namespace math  {
/********************************************************************************/

/** \brief Class for large integers calculus
 *
 * The IntegerTemplate is implemented as a boost variant, which means that it can act like T1, T2, T3 or T4, etc.. type
 * according to the configuration.
 *
 * All the methods are implemented through a boost variant visitor.
 *
 *  Note that we have 2 possible native implementations (NativeInt64 and NativeInt128)
 *  that rely on native types uint64_t and __uint128_t.
 *
 *  For larger integer, a multi-precision LargeInt is used.
 *
 *  From the user point of view, [s]he has just to include this file and use the Integer
 *  class.
 *
 */
template <typename IntegerList>
class IntegerTemplate
{
private:

    /** We define a transformation of the provided integer list in order to get a LargeInt type list. */
    template<typename T>  struct ToLargeInt  {  typedef LargeInt<((T::value+31)/32)> type;  };
    typedef typename boost::mpl::transform<IntegerList, ToLargeInt<boost::mpl::_> >::type transfolist;

public:

    /** We define a boost variant from this type list. */
    typedef typename boost::make_variant_over<transfolist>::type Type;
    /* so, before, Type was defined as a plain boost::variant like this: 
     typedef typename boost::variant<LargeInt<1>, LargeInt<2>, LargeInt<3>, LargeInt<4> > Type;
     I checked: the new make_variant_over has no impact on runtime performance. Probably not surprising, as this code is probably purely compile-time.
     */

    /** Apply a functor with the best template specialization according to the provided kmer size. */
    template <template<size_t> class Functor, typename Parameter>
    static void apply (size_t kmerSize, Parameter params)
    {
        typedef typename boost::mpl::empty<IntegerList>::type empty;

        /** We delegate the execution to the Apply structure, defined with two template specializations
         * that allows recursion. */
        Apply<Functor, Parameter, IntegerList, empty::value>::execute (kmerSize, params);
    }

    /** Constructor. Note that the type (see getType and setType) has to be first initialized
     * otherwise no instance can be created (exception thrown).
     * \param[in] n : value for initialization of the integer.
     */
    IntegerTemplate (int64_t n=0)  {}

    /** Copy constructor. Relies on the copy constructor of boost variant
     * \param[in] t : the object to be used for initialization
     */
    template<typename T>  explicit IntegerTemplate (const T& t) : v (t)  {}

    /** Affectation operator. Relies on the affectation operator of boost variant
     * \param[in] t : object to be copied
     * \return the current object.
     */
    template<typename T>
    IntegerTemplate& operator=(const T& t)
    {
        v = t;
        return *this;
    }

    /** Get the name of the class used by the variant (ie. one of the Ti template class parameters)
     * \return the class name.
     */
    const char*  getName ()         { return boost::apply_visitor (Integer_name(),  *(*this)); }

    /** Get the size of an instance of the class used by the variant  (ie. one of the Ti template class parameters)
     * \return the size of an object (in bits).
     */
    const size_t getSize ()         { return boost::apply_visitor (Integer_size(),  *(*this)); }

    /** Operator +
     * \param[in] a : first operand
     * \param[in] b : second operand
     * \return sum of the two operands.
     */
    inline friend IntegerTemplate   operator+  (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_plus(),  *a, *b);  }

    /** Operator -
     * \param[in] a : first operand
     * \param[in] b : second operand
     * \return substraction of the two operands.
     */
    inline friend IntegerTemplate   operator-  (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_minus(), *a, *b);  }

    /** Operator |
     * \param[in] a : first operand
     * \param[in] b : second operand
     * \return 'or' of the two operands.
     */
    inline friend IntegerTemplate   operator|  (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_or(),    *a, *b);  }

    /** Operator ^
     * \param[in] a : first operand
     * \param[in] b : second operand
     * \return 'xor' of the two operands.
     */
    inline friend IntegerTemplate   operator^  (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_xor(),   *a, *b);  }

    /** Operator &
     * \param[in] a : first operand
     * \param[in] b : second operand
     * \return 'and' of the two operands.
     */
    inline friend IntegerTemplate   operator&  (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_and(),   *a, *b);  }

    /** Operator ~
     * \param[in] a : operand
     * \return negation of the operand
     */
    inline friend IntegerTemplate   operator~  (const IntegerTemplate& a)                    {  Integer_compl v; return  boost::apply_visitor (v, *a);  }

    /** Operator ==
     * \param[in] a : first operand
     * \param[in] b : second operand
     * \return equality of the two operands.
     */
    inline friend bool      operator== (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_equals(), *a, *b);  }

    /** Operator !=
     * \param[in] a : first operand
     * \param[in] b : second operand
     * \return inequality of the two operands.
     */
    inline friend bool      operator!= (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  ! (*a==*b);  }

    /** Operator <
     * \param[in] a : first operand
     * \param[in] b : second operand
     * \return '<' of the two operands.
     */
    inline friend bool      operator<  (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_less(),   *a, *b);  }

    /** Operator <=
     * \param[in] a : first operand
     * \param[in] b : second operand
     * \return '<=' of the two operands.
     */
    inline friend bool      operator<= (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_lesseq(), *a, *b);  }

    /** Operator *
     * \param[in] a : first operand
     * \param[in] c : second operand
     * \return multiplication of the two operands.
     */
    inline friend IntegerTemplate   operator*  (const IntegerTemplate& a, const int&       c)  {  return  boost::apply_visitor (Integer_mult(c), *a);  }

    /** Operator /
     * \param[in] a : first operand
     * \param[in] c : second operand
     * \return division of the two operands.
     */
    inline friend IntegerTemplate   operator/  (const IntegerTemplate& a, const u_int32_t& c)  {  return  boost::apply_visitor (Integer_div(c),  *a);  }

    /** Operator %
     * \param[in] a : first operand
     * \param[in] c : second operand
     * \return modulo of the two operands.
     */
    inline friend u_int32_t operator%  (const IntegerTemplate& a, const u_int32_t& c)  {  return  boost::apply_visitor (Integer_mod(c),  *a);  }

    /** Operator >>
     * \param[in] a : first operand
     * \param[in] c : second operand
     * \return right shift of the two operands.
     */
    inline friend IntegerTemplate   operator>> (const IntegerTemplate& a, const int& c)  {  return  boost::apply_visitor (Integer_shiftLeft(c),   *a);  }

    /** Operator <<
     * \param[in] a : first operand
     * \param[in] c : second operand
     * \return left shift of the two operands.
     */
    inline friend IntegerTemplate   operator<< (const IntegerTemplate& a, const int& c)  {  return  boost::apply_visitor (Integer_shiftRight(c),  *a);  }

    /** Operator +=
     * \param[in] a : first operand
     * \return addition and affectation.
     */
    void operator+= (const IntegerTemplate& a)  {  boost::apply_visitor (Integer_plusaffect(),  *(*this), *a);  }

    /** Operator ^=
     * \param[in] a : first operand
     * \return xor and affectation.
     */
    void operator^= (const IntegerTemplate& a)  {  boost::apply_visitor (Integer_xoraffect(),   *(*this), *a);  }

    /** Operator[] access the ith nucleotide in the given integer. For instance a[4] get the 5th nucleotide of
     * a kmer encoded as an Integer object.
     * \param[in] idx : index of the nucleotide to be retrieved
     * \return the nucleotide value as follow: A=0, C=1, T=2 and G=3
     */
    u_int8_t  operator[]  (size_t idx) const   { return  boost::apply_visitor (Integer_value_at(idx), *(*this)); }

    /** Get the reverse complement of a kmer encoded as an IntegerTemplate object. Note that the kmer size must be known.
     * \param[in] a : kmer value to be reversed-complemented
     * \param[in] sizeKmer : size of the kmer
     * \return the reverse complement kmer as a IntegerTemplate value
     */
    friend IntegerTemplate revcomp (const IntegerTemplate& a,  size_t sizeKmer)  {  return  boost::apply_visitor (Integer_revcomp(sizeKmer),  *a);  }

    /** Get a hash value on 64 bits for a given IntegerTemplate object.
     * \param[in] a : the integer value
     * \param[in] seed : some seed value used for the hash computation.
     * \return the hash value on 64 bits.
     */
    friend u_int64_t hash1        (const IntegerTemplate& a,  u_int64_t seed)  {  return  boost::apply_visitor (Integer_hash1(seed),  *a);          }

    /** Get a hash value on 64 bits for a given IntegerTemplate object.
     * \param[in] a : the integer value
     * \return the hash value on 64 bits.
     */
    friend u_int64_t oahash       (const IntegerTemplate& a)                   {  return  boost::apply_visitor (Integer_oahash(), *a);              }

    /** Get a hash value on 64 bits for a given IntegerTemplate object.
     * Note: although we return 64 bits, only the first 16 bits are set
     * \param[in] a : the integer value
     * \param[in] shift : some value used for the hash computation.
     * \return the hash value on 64 bits.
     */
    friend u_int64_t simplehash16 (const IntegerTemplate& a,  int shift)       {  return  boost::apply_visitor (Integer_simplehash16(shift),  *a);  }

    /** Get the lexicographical minimizer of the k-mer, quickly, using bit tricks. 
     * minimizers containing 'AA' are forbidden. might not return a valid result sometimes.
     * \param[in] nbminimizers: number of minimizers inside this kmer
     * \param[in] &validResult: whether the returned result should be considered valid (if not, compute minimizer with another approach please)
     * \return the reverse complement kmer as a IntegerTemplate value
     */
    friend u_int32_t fastLexiMinimizer (const IntegerTemplate& a, int _nbMinimizers, bool &validResult) {  return  boost::apply_visitor (Integer_fastLexiMinimizer(_nbMinimizers, validResult),  *a);  }



    /** Get an ASCII string representation of a kmer encoded as a IntegerTemplate object
     * \param[in] sizeKmer : size of the kmer
     * \return the ASCII representation of the kmer.
     */
    std::string toString (size_t sizeKmer) const  {  return boost::apply_visitor (Integer_toString(sizeKmer), *(*this)); }

    /** Output stream operator for the IntegerTemplate class
     * \param[in] s : the output stream to be used.
     * \param[in] a : the object to output
     * \return the modified output stream.
     */
    friend std::ostream & operator<<(std::ostream & s, const IntegerTemplate& a)  {  s << *a;  return s;  }

    /** Get the value of the IntegerTemplate object as a U type, U being one of the T1,T2,T3,T4
     * template class parameters. This method can be seen as a converter from the IntegerTemplate class
     * to a specific U type (given as a template parameter of this method).
     * \return  the converted value as a U type.
     */
    template<typename U>
    const U& get ()  const  {  return * boost::get<U>(&v);  }

private:

    struct Integer_name : public boost::static_visitor<const char*>    {
        template<typename T>  const char* operator() (const T& a) const { return a.getName();  }};

    struct Integer_size : public boost::static_visitor<const size_t>    {
        template<typename T>  const size_t operator() (const T& a) const  { return a.getSize();  }};

    struct Integer_plus : public boost::static_visitor<IntegerTemplate>    {
        template<typename T>              IntegerTemplate operator() (const T& a, const T& b) const  { return IntegerTemplate(a + b);  }
        template<typename T, typename U>  IntegerTemplate operator() (const T& a, const U& b) const  { return IntegerTemplate();       }
    };

    struct Integer_minus : public boost::static_visitor<IntegerTemplate>    {
        template<typename T>              IntegerTemplate operator() (const T& a, const T& b) const  { return IntegerTemplate(a - b);  }
        template<typename T, typename U>  IntegerTemplate operator() (const T& a, const U& b) const  { return IntegerTemplate();  }
    };

    struct Integer_or : public boost::static_visitor<IntegerTemplate>    {
        template<typename T>              IntegerTemplate operator() (const T& a, const T& b) const  { return IntegerTemplate(a | b);  }
        template<typename T, typename U>  IntegerTemplate operator() (const T& a, const U& b) const  { return IntegerTemplate();  }
    };

    struct Integer_xor : public boost::static_visitor<IntegerTemplate>    {
        template<typename T>              IntegerTemplate operator() (const T& a, const T& b) const  { return IntegerTemplate(a ^ b);  }
        template<typename T, typename U>  IntegerTemplate operator() (const T& a, const U& b) const  { return IntegerTemplate();  }
    };

    struct Integer_and : public boost::static_visitor<IntegerTemplate>    {
        template<typename T>              IntegerTemplate operator() (const T& a, const T& b) const  { return IntegerTemplate(a & b);  }
        template<typename T, typename U>  IntegerTemplate operator() (const T& a, const U& b) const  { return IntegerTemplate();  }
    };

    struct Integer_less : public boost::static_visitor<bool>    {
        template<typename T>              bool operator() (const T& a, const T& b) const  { return a < b;  }
        template<typename T, typename U>  bool operator() (const T& a, const U& b) const  { return false;  }
    };

    struct Integer_lesseq : public boost::static_visitor<bool>    {
        template<typename T>              bool operator() (const T& a, const T& b) const  { return a <= b;  }
        template<typename T, typename U>  bool operator() (const T& a, const U& b) const  { return false;   }
    };

    struct Integer_equals : public boost::static_visitor<bool>    {
        template<typename T>              bool operator() (const T& a, const T& b) const  { return a == b;  }
        template<typename T, typename U>  bool operator() (const T& a, const U& b) const  { return false;   }
    };

    struct Integer_plusaffect : public boost::static_visitor<>    {
        template<typename T>              void operator() ( T& a, const T& b) const  { a += b;  }
        template<typename T, typename U>  void operator() ( T& a, const U& b) const  {   }
    };

    struct Integer_xoraffect : public boost::static_visitor<>    {
        template<typename T>              void operator() ( T& a, const T& b) const  { a ^= b;  }
        template<typename T, typename U>  void operator() ( T& a, const U& b) const  {   }
    };

    struct Integer_compl : public boost::static_visitor<IntegerTemplate>    {
        template<typename T>  IntegerTemplate operator() (const T& a)  { return IntegerTemplate(~a);  }};

    template<typename Result, typename Arg>
    struct Visitor : public boost::static_visitor<Result>
    {
        Visitor (Arg a=Arg()) : arg(a) {}
        Arg arg;
    };
	
    template<typename Result, typename Arg1, typename Arg2>
    struct Visitor2Args : public boost::static_visitor<Result>
    {
        Visitor2Args (Arg1 a1=Arg1(), Arg2 a2=Arg2()) : arg1(a1), arg2(a2) {}
        Arg1 arg1;
        Arg2 arg2;
    };

    struct Integer_mult : public Visitor<IntegerTemplate,const int>    {
        Integer_mult (const int& c) : Visitor<IntegerTemplate,const int>(c) {}
        template<typename T>  IntegerTemplate operator() (const T& a) const  { return IntegerTemplate(a*this->arg);  }};

    struct Integer_div : public Visitor<IntegerTemplate,const u_int32_t>    {
        Integer_div (const u_int32_t& c) : Visitor<IntegerTemplate,const u_int32_t>(c) {}
        template<typename T>  IntegerTemplate operator() (const T& a) const  { return IntegerTemplate(a/this->arg);  }};

    struct Integer_mod : public Visitor<u_int32_t,const u_int32_t>    {
        Integer_mod (const u_int32_t& c) : Visitor<u_int32_t,const u_int32_t>(c) {}
        template<typename T>  u_int32_t operator() (const T& a) const  { return (a%this->arg);  }};

    struct Integer_shiftLeft : public Visitor<IntegerTemplate,const int>    {
        Integer_shiftLeft (const int& c) : Visitor<IntegerTemplate,const int>(c) {}
        template<typename T>  IntegerTemplate operator() (const T& a) const  { return IntegerTemplate (a >> this->arg);  }};

    struct Integer_shiftRight : public Visitor<IntegerTemplate,const int>    {
        Integer_shiftRight (const int& c) : Visitor<IntegerTemplate,const int>(c) {}
        template<typename T>  IntegerTemplate operator() (const T& a) const  { return IntegerTemplate (a << this->arg);  }};

    struct Integer_revcomp : public Visitor<IntegerTemplate,size_t>    {
        Integer_revcomp (const size_t& c) : Visitor<IntegerTemplate,size_t>(c) {}
        template<typename T>  IntegerTemplate operator() (const T& a) const  { return IntegerTemplate (revcomp(a,this->arg));  }};

    struct Integer_hash1 : public Visitor<u_int64_t,u_int64_t>    {
        Integer_hash1 (const u_int64_t& c) : Visitor<u_int64_t,u_int64_t>(c) {}
        template<typename T>  u_int64_t operator() (const T& a) const  { return (hash1(a,this->arg));  }};

    struct Integer_oahash : public boost::static_visitor<u_int64_t>    {
        template<typename T>  u_int64_t operator() (const T& a) const  { return (oahash(a));  }};

    struct Integer_simplehash16 : public Visitor<u_int64_t,int>    {
        Integer_simplehash16 (const int& c) : Visitor<u_int64_t,int>(c) {}
        template<typename T>  u_int64_t operator() (const T& a) const  { return (simplehash16(a,this->arg));  }};

    struct Integer_fastLexiMinimizer : public Visitor2Args<u_int32_t, int, bool&>    {
        Integer_fastLexiMinimizer (int b, bool& c) : Visitor2Args<u_int32_t,int, bool&>(b,c) {}
        template<typename T>  u_int32_t operator() (const T& a) const  { return (fastLexiMinimizer (a,this->arg1, this->arg2));  }};
	
    struct Integer_value_at : public Visitor<u_int8_t,size_t>   {
        Integer_value_at (size_t idx) : Visitor<u_int8_t,size_t>(idx) {}
        template<typename T>  u_int8_t operator() (const T& a) const { return a[this->arg];  }};

    struct Integer_toString : public Visitor<std::string,size_t>   {
        Integer_toString (size_t c) : Visitor<std::string,size_t>(c) {}
        template<typename T>  std::string operator() (const T& a) const  { return a.toString(this->arg);  }};

private:

    /** We instantiate the boost variant. */
    Type v;

          Type& operator *()       { return v; }
    const Type& operator *() const { return v; }


    /** Now, we define the Apply structure that allows to find the correct implementation of LargeInt
     * according to the given kmerSize (at runtime). */

    /** My initial guess didn't work, although it should have (pb with boost::mpl ?)... I went back on an implementation
     * similar to http://www.developpez.net/forums/d1193120/c-cpp/cpp/bibliotheques/boost/vecteur-vide-boost-mpl */

    /** Template definition. */
    template<template<size_t> class Functor, class Parameter, class T, bool empty> struct Apply  {};

    /** Template specialization for a vector of types. */
    template<template<size_t> class Functor, class Parameter, class T>
    struct Apply<Functor, Parameter, T, false>
    {
        static void execute (size_t kmerSize, Parameter params)
        {
            /** Shortcut : we get the current kmer size threshold. */
            static const size_t K = boost::mpl::front<T>::type::value;

            /** We check whether the kmerSize parameter falls into the current K threshold.
             * If yes, we run the functor and leave. */
            if (kmerSize < K)  { Functor<K>() (params);  return; }

            typedef typename boost::mpl::pop_front<T>::type tail;
            typedef typename boost::mpl::empty<tail>::type  empty;

            /** The current K threshold doesn't work, try the next kmer value. */
            Apply<Functor, Parameter, tail, empty::value>::execute (kmerSize, params);
        }
    };

    /** Template specialization for an empty type. */
    template<template<size_t> class Functor, class Parameter, class T>
    struct Apply<Functor, Parameter, T, true>
    {
        static void execute (size_t kmerSize, Parameter params)
        {
            throw system::Exception ("Failure because of unhandled kmer size %d", kmerSize);
        }
    };
};


// this is just for benchmarking. It helped me see that indeed, instantiating a class that constains a boost::variant takes significant run-time :(
template <typename IntegerList>
class IntegerTemplateDummy
{
public:

    typedef typename boost::variant<LargeInt<1>, LargeInt<2>, LargeInt<3>, LargeInt<4> > Type;

private:
    /** We instantiate the boost variant. */
    Type v;
};

/********************************************************************************/

/** We define a mpl::vector holding the int_<K> types, one per kmer size value chosen by the user. */
typedef boost::mpl::vector<KSIZE_LIST_TYPE>::type  IntegerList;

/** We specialize IntegerTemplate class based on the list of kmer size values chosen by the user. */
typedef  IntegerTemplate <IntegerList> Integer;
typedef  IntegerTemplateDummy<IntegerList> IntegerDummy;

/********************************************************************************/
}}}};
/********************************************************************************/

/* would enable an Integer to be hashed in a unordered_map */
/* well, for some reason this didn't work with unordered_map, i still got the compile error:
 /usr/include/c++/6.1.1/bits/hashtable_policy.h:85:34: error: no match for call to ‘(const std::hash<const gatb::core::tools::math::IntegerTemplate<boost::mpl::vector1<mpl_::int_<32> > > >) (const gatb::core::tools::math::IntegerTemplate<boost::mpl::vector1<mpl_::int_<32> > >&)’ */
/*namespace std {
  template <typename IntegerList>
  struct hash< gatb::core::tools::math::IntegerTemplate<IntegerList>>
  {
    std::size_t operator()(const gatb::core::tools::math::IntegerTemplate<IntegerList>& k) const
    {
        return oahash(k);
    }
  };
};
*/

#endif /* _GATB_CORE_TOOLS_MATH_INTEGER_HPP_ */
