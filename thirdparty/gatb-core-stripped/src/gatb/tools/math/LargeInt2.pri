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

/** \file LargeInt<2>.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Integer class relying on native u_int64_t type
 */

/********************************************************************************/
#ifdef __SIZEOF_INT128__
/********************************************************************************/

template<>  class LargeInt<2> 
{
public:

#ifdef USE_LARGEINT_CONSTRUCTOR 
    /** Constructor.
     * \param[in] c : initial value of the large integer. */
    LargeInt<2> (const __uint128_t& c=0)  {  value = c;  }
#endif

     u_int64_t getVal () const  { return value; }
     inline void setVal (const u_int64_t &c) { value = c; }
     inline void setVal (const LargeInt<2>& c) { value = c.value; }

    static const char* getName ()  { return "LargeInt<2>"; }

    static const size_t getSize ()  { return 8*sizeof(__uint128_t); }

    LargeInt<2> operator+  (const LargeInt<2>& other)   const   {  LargeInt<2> res; res.value = value + other.value; return res; }
    LargeInt<2> operator+  (const u_int64_t& other)     const   {  LargeInt<2> res; res.value = value + other; return res; }
    LargeInt<2> operator-  (const LargeInt<2>& other)   const   {  LargeInt<2> res; res.value = value - other.value; return res; }
    LargeInt<2> operator-  (const u_int64_t& other)     const   {  LargeInt<2> res; res.value = value - other; return res; }
    LargeInt<2> operator|  (const LargeInt<2>& other)   const   {   LargeInt<2> res; res.value = value | other.value; return res; }
    LargeInt<2> operator*  (const int& coeff)           const   {   LargeInt<2> res; res.value = value * coeff;       return res; }
    LargeInt<2> operator/  (const u_int32_t& divisor)   const   {   LargeInt<2> res; res.value = value / divisor;     return res; }
    u_int32_t   operator%  (const u_int32_t& divisor)   const   {   return value % divisor;  }
    LargeInt<2> operator^  (const LargeInt<2>& other)   const   {   LargeInt<2> res; res.value = value ^ other.value; return res; }
    LargeInt<2> operator&  (const LargeInt<2>& other)   const   {   LargeInt<2> res; res.value = value & other.value; return res; }
    LargeInt<2> operator&  (const char& other)          const   {   LargeInt<2> res; res.value = value & other;       return res; }
    LargeInt<2> operator~  ()                           const   {   LargeInt<2> res; res.value = ~value;              return res; }
    LargeInt<2> operator<< (const int& coeff)           const   {   LargeInt<2> res; res.value = value << coeff;      return res; }
    LargeInt<2> operator>> (const int& coeff)           const   {   LargeInt<2> res; res.value = value >> coeff;      return res; }
    /*    LargeInt<2> operator=  (const u_int64_t& other)     const   {   LargeInt<2> res; res.value = other; return res; }
    LargeInt<2> operator=  (const LargeInt<2>& other)   const   {   LargeInt<2> res; res.value = other.value; return res; }*/ // that code is wrong. see LargeInt1.pri for a fix for operator=(const LargeInt&), haven't figured out for operator=(uint64), decided to use setVal instead

    bool         operator!= (const LargeInt<2>& c)         const   {  return value != c.value;     }
    bool         operator!= (const u_int64_t  & c)         const   {  return value != c;           }
    bool         operator== (const LargeInt<2>& c)         const   {  return value == c.value;     }
    bool         operator== (const u_int64_t  & c)         const   {  return value == c;           }
    bool         operator<  (const LargeInt<2>& c)         const   {  return value < c.value;      }
    bool         operator<= (const LargeInt<2>& c)         const   {  return value <= c.value;     }

    LargeInt<2>& operator+=  (const LargeInt<2>& other)    {  value += other.value; return *this; }
    LargeInt<2>& operator+=  (const u_int64_t& other)      {  value += other; return *this; }
    LargeInt<2>& operator^=  (const LargeInt<2>& other)    {  value ^= other.value; return *this; }

    LargeInt<2>& operator<<=  (const int& coeff)  { value <<= coeff; return *this; } 
    LargeInt<2>& operator>>=  (const int& coeff)  { value >>= coeff; return *this; }

    u_int8_t  operator[]  (size_t idx) const   {  return (value >> (2*idx)) & 3; }

    /** Output stream overload. NOTE: for easier process, dump the value in hexadecimal.
     * \param[in] os : the output stream
     * \param[in] in : the integer value to be output.
     * \return the output stream.
     */
    friend std::ostream & operator<<(std::ostream & os, const LargeInt<2> & in)
    {
        __uint128_t x = in.value;

        u_int64_t high_nucl = (u_int64_t) (x>>64);
        u_int64_t low_nucl  = (u_int64_t)(x&((((__uint128_t)1)<<64)-1));

        if (high_nucl == 0) {   os << std::hex <<                     low_nucl << std::dec;  }
        else                {   os << std::hex << high_nucl << "." << low_nucl << std::dec;  }
        return os;
    }

    /********************************************************************************/
    
    /** Print corresponding kmer in ASCII
     * \param[sizeKmer] in : kmer size (def=64).
     */
    inline void printASCII ( size_t sizeKmer = 64)
    {
        int i;
        u_int64_t temp = value;
        
        
        char seq[65];
        char bin2NT[4] = {'A','C','T','G'};
        
        for (i=sizeKmer-1; i>=0; i--)
        {
            seq[i] = bin2NT[ temp&3 ];
            temp = temp>>2;
        }
        seq[sizeKmer]='\0';
        
        std::cout << seq << std::endl;
    }
    
    /********************************************************************************/
    /** Print corresponding kmer in ASCII
     * \param[sizeKmer] in : kmer size (def=32).
     */
    std::string toString (size_t sizeKmer) const
    {
        char seq[65];
        char bin2NT[4] = {'A','C','T','G'};

        for (size_t i=0; i<sizeKmer; i++)  {  seq[sizeKmer-i-1] = bin2NT [(*this)[i]];  }
        seq[sizeKmer]='\0';
        return seq;
    }

    /********************************************************************************/
    template<typename Map>
    static LargeInt<2> polynom (const char* data, size_t size, Map fct)
    {
        LargeInt<2> res;
        res.setVal(0);
        for (size_t i=0; i<size; ++i)  {  res.value = res.value * 4 + fct(data[i]);  }
        return res;
    }

    __uint128_t get_128() const
    {
        return value;
    }

    const uint64_t* get_data() const
    {
        return reinterpret_cast<const uint64_t*>(&value);
    }

private:
    __uint128_t value;
    friend LargeInt<2> revcomp (const LargeInt<2>& i,   size_t sizeKmer);
    friend u_int64_t    hash1    (const LargeInt<2>& key, u_int64_t  seed);
    friend u_int64_t    hash2    (const LargeInt<2>& key, u_int64_t  seed);
    friend u_int64_t    oahash  (const LargeInt<2>& key);
    friend u_int64_t    simplehash16    (const LargeInt<2>& key, int  shift);
    template<typename m_T> friend void fastLexiMinimizer (const LargeInt<2>& x, const unsigned int _nbMinimizers, \
                             const unsigned int m, m_T &minimizer, size_t &position, bool &validResult);
    friend void justSweepForAA(const LargeInt<2>& x, const unsigned int _nbMinimizers, unsigned int &dummy);

};

/********************************************************************************/
inline LargeInt<2> revcomp (const LargeInt<2>& in, size_t sizeKmer)
{
    //                  ---64bits--   ---64bits--
    // original kmer: [__high_nucl__|__low_nucl___]
    //
    // ex:            [         AC  | .......TG   ]
    //
    //revcomp:        [         CA  | .......GT   ]
    //                 \_low_nucl__/\high_nucl/

    const __uint128_t& x = in.value;

    u_int64_t high_nucl = (u_int64_t)(x>>64);
    int nb_high_nucl = sizeKmer>32?sizeKmer - 32:0;

    u_int64_t revcomp_high_nucl = NativeInt64::revcomp64 (high_nucl, nb_high_nucl);

    if (sizeKmer<=32) revcomp_high_nucl = 0; // srsly dunno why this is needed. gcc bug? u_int64_t x ---> (x>>64) != 0

    u_int64_t low_nucl = (u_int64_t)(x&((((__uint128_t)1)<<64)-1));
    int nb_low_nucl = sizeKmer>32?32:sizeKmer;

    u_int64_t revcomp_low_nucl = NativeInt64::revcomp64 (low_nucl, nb_low_nucl);

    LargeInt<2> res;
    res.setVal(revcomp_low_nucl);
    res <<= 2* nb_high_nucl;
    res += revcomp_high_nucl;
    return res;
}

/********************************************************************************/
inline u_int64_t hash1 (const LargeInt<2>& item, u_int64_t seed=0)
{
    const __uint128_t& elem = item.value;

    return NativeInt64::hash64 ((u_int64_t)(elem>>64),seed) ^
           NativeInt64::hash64 ((u_int64_t)(elem&((((__uint128_t)1)<<64)-1)),seed);
}

inline u_int64_t hash2 (const LargeInt<2>& item, u_int64_t seed=0)
{
    const __uint128_t& elem = item.value;

        // from inline uint64_t twang_mix64(uint64_t key) taken from https://github.com/facebook/folly/blob/master/folly/Hash.h
        u_int64_t key = (u_int64_t)(elem>>64);
        key = (~key) + (key << 21);  // key *= (1 << 21) - 1; key -= 1;
        key = key ^ (key >> 24);
        key = key + (key << 3) + (key << 8);  // key *= 1 + (1 << 3) + (1 << 8)
        key = key ^ (key >> 14);
        key = key + (key << 2) + (key << 4);  // key *= 1 + (1 << 2) + (1 << 4)
        key = key ^ (key >> 28);
        key = key + (key << 31);  // key *= 1 + (1 << 31)
        u_int64_t res = key;

        key = (u_int64_t)(elem&((((__uint128_t)1)<<64)-1));
        key = (~key) + (key << 21);  // key *= (1 << 21) - 1; key -= 1;
        key = key ^ (key >> 24);
        key = key + (key << 3) + (key << 8);  // key *= 1 + (1 << 3) + (1 << 8)
        key = key ^ (key >> 14);
        key = key + (key << 2) + (key << 4);  // key *= 1 + (1 << 2) + (1 << 4)
        key = key ^ (key >> 28);
        key = key + (key << 31);  // key *= 1 + (1 << 31)

    return res ^ key;
}


/********************************************************************************/
inline u_int64_t oahash (const LargeInt<2>& item)
{
    const __uint128_t& elem = item.value;

    return NativeInt64::oahash64 ((u_int64_t)(elem>>64)) ^
           NativeInt64::oahash64 ((u_int64_t)(elem&((((__uint128_t)1)<<64)-1)));
}



/********************************************************************************/
inline u_int64_t simplehash16 (const LargeInt<2>& key, int  shift)
{
    return NativeInt64::simplehash16_64 ((u_int64_t)key.value, shift);
}

/********************************************************************************/
template<typename minimizer_type> void fastLexiMinimizer (const LargeInt<2>& x, const unsigned int _nbMinimizers, \
                             const unsigned int m, minimizer_type &minimizer, size_t &position, bool &validResult)
{
    const minimizer_type default_minimizer = ~0 & ((1 << (2*m)) - 1); 
    minimizer = default_minimizer; 
    validResult = false;

    fastLexiMinimizerChunk<__uint128_t,minimizer_type>(x.value, _nbMinimizers, m, 0,  minimizer, position, 0, validResult);

    validResult = validResult && (minimizer != default_minimizer) /* might happen that AA was found but resulted in forbidden minimizers */;
}

/* debug function, for profiling only; counts the AA's in a kmer */
inline void justSweepForAA(const LargeInt<2>& x, const unsigned int _nbMinimizers, unsigned int &dummy) 
{
        __uint128_t val = x.value;

        const int it = std::min((unsigned int)sizeof(__uint128_t)*4, _nbMinimizers); 
        int j = 0;
        while (j < it)
        {
            if ((val & 15) == 0) // val starts with AA
                dummy++;

            val >>= 2;
        }
}


/********************************************************************************/
#endif
/********************************************************************************/
