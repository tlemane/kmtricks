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
/** \file LargeInt<1>.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Integer class relying on native u_int64_t type
 */


template<>  class LargeInt<1> 
{
public:

#ifdef USE_LARGEINT_CONSTRUCTOR
    /** Constructor.
     * \param[in] c : initial value of the large integer. */
    LargeInt<1> (const u_int64_t& c=0)  {  value = c;  }
#endif

// due to a bug in ograph, reintroducing constructors provided a quick fix. here they are. I think the problem was caused because i added the operator=()'s 
#if 0
    LargeInt<1> ()  {  value = 0; }
    LargeInt<1> (const u_int64_t& c)  {  value = c;  }
    LargeInt<1> (const LargeInt<1>& other)  {  value = other.value;  }
#endif

     u_int64_t getVal () const   { return value; }
     inline void setVal (u_int64_t val) { value = val; }
     inline void setVal (const LargeInt<1>& other) { value = other.value; }

    static const char* getName ()  { return "LargeInt<1>"; }

    static const size_t getSize ()  { return 8*sizeof(u_int64_t); }

    /** Returns lower 64 bits */
    u_int64_t toInt () const  {  return value;  }

    LargeInt<1> operator+  (const LargeInt<1>& other)   const   {  LargeInt<1> res; res.value = value + other.value; return res; }
    LargeInt<1> operator+  (const u_int64_t& other)     const   {  LargeInt<1> res; res.value = value + other; return res; }
    LargeInt<1> operator-  (const LargeInt<1>& other)   const   {  LargeInt<1> res; res.value = value - other.value; return res; }
    LargeInt<1> operator-  (const u_int64_t& other)     const   {  LargeInt<1> res; res.value = value - other; return res; }
    LargeInt<1> operator|  (const LargeInt<1>& other)   const   {   LargeInt<1> res; res.value = value | other.value; return res; }
    LargeInt<1> operator*  (const int& coeff)           const   {   LargeInt<1> res; res.value = value * coeff;       return res; }
    LargeInt<1> operator/  (const u_int32_t& divisor)   const   {   LargeInt<1> res; res.value = value / divisor;     return res; }
    u_int32_t   operator%  (const u_int32_t& divisor)   const   {   return value % divisor;  }
    LargeInt<1> operator^  (const LargeInt<1>& other)   const   {   LargeInt<1> res; res.value = value ^ other.value; return res; }
    LargeInt<1> operator&  (const LargeInt<1>& other)   const   {   LargeInt<1> res; res.value = value & other.value; return res; }
    LargeInt<1> operator&  (const char& other)          const   {   LargeInt<1> res; res.value = value & other;       return res; }
    LargeInt<1> operator~  ()                           const   {   LargeInt<1> res; res.value = ~value;              return res; }
    LargeInt<1> operator<< (const int& coeff)           const   {   LargeInt<1> res; res.value = value << coeff;      return res; }
    LargeInt<1> operator>> (const int& coeff)           const   {   LargeInt<1> res; res.value = value >> coeff;      return res; }
#if 0  // last time i defined those copy operators, i ran into troubles, so let's not do that and let's use setVal instead!
    LargeInt<1> operator=  (const uint64_t& c)          const   {   res.value = c;  return res; }
    LargeInt<1> operator=  (const LargeInt<1>& other)   const   {   if (this == &other) return *this; res.value = other.value;  return *this; }
#endif
    bool        operator!= (const LargeInt<1>& c)       const   {   return value != c.value; }
    bool        operator!= (const u_int64_t& c)         const   {   return value != c; }
    bool        operator== (const LargeInt<1>& c)       const   {   return value == c.value; }
    bool        operator== (const u_int64_t& c)         const   {   return value == c; }
    bool        operator<  (const LargeInt<1>& c)       const   {   return value < c.value; }
    bool        operator<= (const LargeInt<1>& c)       const   {   return value <= c.value; }

    LargeInt<1>& operator+=  (const LargeInt<1>& other)    {  value += other.value; return *this; }
    LargeInt<1>& operator^=  (const LargeInt<1>& other)    {  value ^= other.value; return *this; }

    LargeInt<1>& operator<<=  (const int& coeff)  { value <<= coeff; return *this; } 
    LargeInt<1>& operator>>=  (const int& coeff)  { value >>= coeff; return *this; }

    u_int8_t  operator[]  (size_t idx) const   {  return (value >> (2*idx)) & 3; }

    /********************************************************************************/
    friend std::ostream & operator<<(std::ostream & s, const LargeInt<1> & l)
    {
        s << std::hex << l.value << std::dec;  return s;
    }
    /********************************************************************************/
    /** Print corresponding kmer in ASCII
     * \param[sizeKmer] in : kmer size (def=32).
     */
    inline void printASCII ( size_t sizeKmer = 32)
    {
        int i;
        u_int64_t temp = value;

        
        char seq[33];
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
        int i;
        u_int64_t temp = value;

        char seq[33];
        char bin2NT[4] = {'A','C','T','G'};

        for (i=sizeKmer-1; i>=0; i--)
        {
            seq[i] = bin2NT[ temp&3 ];
            temp = temp>>2;
        }
        seq[sizeKmer]='\0';
        return seq;
    }

    /********************************************************************************/
    inline static u_int64_t revcomp64 (const u_int64_t& x, size_t sizeKmer)
    {
        u_int64_t res = x;

        // OLD VERSION (with lookup table)
        // unsigned char* kmerrev  = (unsigned char *) (&(res));
        // unsigned char* kmer     = (unsigned char *) (&(x));
        // for (size_t i=0; i<8; ++i)  {  kmerrev[8-1-i] = revcomp_4NT [kmer[i]];  }

        res = ((res>> 2 & 0x3333333333333333) | (res & 0x3333333333333333) <<  2);
        res = ((res>> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) <<  4);
        res = ((res>> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) <<  8);
        res = ((res>>16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
        res = ((res>>32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
        res = res ^ 0xAAAAAAAAAAAAAAAA;
        
        return (res >> (2*(32-sizeKmer))) ;
    }

    /********************************************************************************/
    inline static u_int64_t hash64 (u_int64_t key, u_int64_t seed)
    {
        u_int64_t hash = seed;
        hash ^= (hash <<  7) ^  key * (hash >> 3) ^ (~((hash << 11) + (key ^ (hash >> 5))));
        hash = (~hash) + (hash << 21); // hash = (hash << 21) - hash - 1;
        hash = hash ^ (hash >> 24);
        hash = (hash + (hash << 3)) + (hash << 8); // hash * 265
        hash = hash ^ (hash >> 14);
        hash = (hash + (hash << 2)) + (hash << 4); // hash * 21
        hash = hash ^ (hash >> 28);
        hash = hash + (hash << 31);
		
		return hash;
    }

    /********************************************************************************/
    inline static u_int64_t oahash64 (u_int64_t elem)
    {
        u_int64_t code = elem;
        code = code ^ (code >> 14); //supp
        code = (~code) + (code << 18);
        code = code ^ (code >> 31);
        code = code * 21;
        code = code ^ (code >> 11);
        code = code + (code << 6);
        code = code ^ (code >> 22);
        return code;
    }

    /********************************************************************************/
    /** computes a simple, naive hash using only 16 bits from input key
     * \param[shift] in : selects which of the input byte will be used for hash computation
     */
    inline static  u_int64_t    simplehash16_64   (u_int64_t key, int  shift)
    {
        u_int64_t input = key >> shift;
        u_int64_t res = random_values[input & 255]   ;
        
        input = input  >> 8;
        res  ^= random_values[input & 255] ;

        
        res  ^= random_values[key & 255] ;//also always add 8 first bits
//		
//		input = input  >> 8;
//        res  ^= random_values[input & 255] ;
//		
//		
//		input = input  >> 8;
//        res  ^= random_values[input & 255] ;
		
		
        return res;
        //could be improved by xor'ing result of multiple bytes
    }

    /********************************************************************************/
    template<typename Map>
    static LargeInt<1> polynom (const char* data, size_t size, Map fct)
    {
        LargeInt<1> res;
        res.value = 0;
        for (size_t i=0; i<size; ++i)  {  res.value = 4 * res.value + fct(data[i]);  }
        return res;
    }

    u_int64_t value; // not ArrayData anymore

    const uint64_t* get_data() const
    {
        return &value;
    }

private:

    friend LargeInt<1> revcomp (const LargeInt<1>& i,   size_t sizeKmer);
    friend u_int64_t    hash1    (const LargeInt<1>& key, u_int64_t  seed);
    friend u_int64_t    hash2    (const LargeInt<1>& key, u_int64_t  seed);
    friend u_int64_t    oahash  (const LargeInt<1>& key);
    friend u_int64_t    simplehash16    (const LargeInt<1>& key, int  shift);
    friend void fastLexiMinimizer (const LargeInt<1>& x, const unsigned int _nbMinimizers, const unsigned int m,  u_int32_t &minimizer, size_t &position, bool &validResult);

};

/********************************************************************************/
inline LargeInt<1> revcomp (const LargeInt<1>& x, size_t sizeKmer)
{
    LargeInt<1> res;
    res.value = LargeInt<1>::revcomp64 (x.value, sizeKmer);
    return res;
}

/********************************************************************************/
inline u_int64_t hash1 (const LargeInt<1>& key, u_int64_t seed=0)
{

    return LargeInt<1>::hash64 (key.value, seed);
}

inline u_int64_t hash2 (const LargeInt<1>& x, u_int64_t seed=0)
{
// from inline uint64_t twang_mix64(uint64_t key) taken from https://github.com/facebook/folly/blob/master/folly/Hash.h
  u_int64_t key = x.value;
  key = (~key) + (key << 21);  // key *= (1 << 21) - 1; key -= 1;
  key = key ^ (key >> 24);
  key = key + (key << 3) + (key << 8);  // key *= 1 + (1 << 3) + (1 << 8)
  key = key ^ (key >> 14);
  key = key + (key << 2) + (key << 4);  // key *= 1 + (1 << 2) + (1 << 4)
  key = key ^ (key >> 28);
  key = key + (key << 31);  // key *= 1 + (1 << 31)
  return key;
}

/********************************************************************************/
inline u_int64_t oahash (const LargeInt<1>& key)
{
    return LargeInt<1>::oahash64 (key.value);
}

/********************************************************************************/
inline u_int64_t simplehash16 (const LargeInt<1>& key, int  shift)
{
    return LargeInt<1>::simplehash16_64 (key.value, shift);
}

inline void fastLexiMinimizer (const LargeInt<1>& x, const unsigned int _nbMinimizers, const unsigned int m,  u_int32_t &minimizer, size_t &position, bool &validResult) 
{
    if (m > sizeof(u_int32_t)*4) {std::cout << "wrong minimizer size for fastLeximinimizer :" << m; exit(1);}

    const u_int32_t default_minimizer = ~0 & ((1 << (2*m)) - 1); 
    minimizer = default_minimizer; 

    validResult = false;
    bool AA_found = false;

    u_int64_t val = x.value;
    const u_int32_t high_bits = 0;
    const int position_offset = 0;

    fastLexiMinimizerChunk<u_int64_t,u_int32_t>(val, _nbMinimizers, m, high_bits, minimizer, position, position_offset, AA_found);

    validResult = AA_found && (minimizer != default_minimizer) /* might happen that AA was found but resulted in forbidden minimizers */;
}

/* debug function, for profiling only; counts the AA's in a kmer */
inline void justSweepForAA(const LargeInt<1>& x, const unsigned int _nbMinimizers, unsigned int &dummy) 
{
    u_int64_t val = x.value;

    const int it = std::min((unsigned int)sizeof(u_int64_t)*4,(unsigned int) _nbMinimizers); 
    int j = 0;
    while (j < it)
    {
        if ((val & 15) == 0) // val starts with AA
            dummy++;

        val >>= 2;
        j++;
    }
}

