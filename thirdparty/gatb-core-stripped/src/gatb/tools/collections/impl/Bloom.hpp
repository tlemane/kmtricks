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

/** \file Bloom.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Bloom implementation
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BLOOM_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BLOOM_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Container.hpp>
#include <gatb/tools/collections/api/Bag.hpp>
#include <gatb/tools/math/LargeInt.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/system/api/types.hpp>
#include <gatb/tools/misc/api/Enums.hpp>
#include <bitset>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

extern u_int8_t bit_mask [];

/********************************************************************************/
/** \brief Define a set of hash functions
 *
 * This class defines a set of hash functions for a given (template) type.
 *
 * One can get a hash code for a given [function,item] through the operator()
 *
 * Note: this class is mainly used by Bloom filters implementations. It is not
 * primarly targeted for end users but they could nevertheless use it.
 */
template <typename Item> class HashFunctors
{
public:

    /** Constructor.
     * \param[in] nbFct : number of hash functions to be used
     * \param[in] seed : some initialization code for defining the hash functions. */
    HashFunctors (size_t nbFct, u_int64_t seed=0) : _nbFct(nbFct), user_seed(seed)
    {
        generate_hash_seed ();
    }

    /** Get a hash code for a hash function and a given item.
     * \param[in] key : item for which we want a hash code
     * \param[in] idx : index of the hash function to be used
     * \return the hash code for the item. */
    u_int64_t operator ()  (const Item& key, size_t idx)  {  return hash1 (key, seed_tab[idx]);  }

private:

    /* */
    void generate_hash_seed ()
    {
        static const u_int64_t rbase[NSEEDSBLOOM] =
        {
            0xAAAAAAAA55555555ULL,  0x33333333CCCCCCCCULL,  0x6666666699999999ULL,  0xB5B5B5B54B4B4B4BULL,
            0xAA55AA5555335533ULL,  0x33CC33CCCC66CC66ULL,  0x6699669999B599B5ULL,  0xB54BB54B4BAA4BAAULL,
            0xAA33AA3355CC55CCULL,  0x33663366CC99CC99ULL
        };

        for (size_t i=0; i<NSEEDSBLOOM; ++i)  {  seed_tab[i] = rbase[i];  }
        for (size_t i=0; i<NSEEDSBLOOM; ++i)  {  seed_tab[i] = seed_tab[i] * seed_tab[(i+3) % NSEEDSBLOOM] + user_seed ;  }
    }

    size_t _nbFct;

    static const size_t NSEEDSBLOOM = 10;
    u_int64_t seed_tab[NSEEDSBLOOM];
    u_int64_t user_seed;
};

/********************************************************************************/

/** \brief Bloom interface
 *
 * We define a Bloom filter has a container (something that tells whether an item is here or not) and
 * a bag (something we can put items into).
 *
 * This interface has a template argument corresponding to the type of items to be inserted in the filter.
 * Note that implementations should provide hash functions supporting the template type.
 *
 * As expected, there is no Iterable interface here because it is not possible to enumerate the items
 * inserted in a Bloom filter only with the Bloom filter information.
 */
template <typename Item> class IBloom : public Container<Item>, public Bag<Item>, public system::SmartPointer
{
public:

    /** Destructor. */
    virtual ~IBloom() {}

    /** Get the raw bit set of the Bloom filter.
     * \return Bloom filter's bit set. */
    virtual u_int8_t*& getArray    () = 0;

    /** Get the size of the Bloom filter (in bytes).
     * \return the size of the bit set of the Bloom filter */
    virtual u_int64_t  getSize     () = 0;

    /** Get the size of the Bloom filter (in bits).
     * \return the size of the bit set of the Bloom filter */
    virtual u_int64_t  getBitSize  () = 0;

    /** Get the number of hash functions used for the Bloom filter.
     * \return the number of hash functions. */
    virtual size_t     getNbHash   () const = 0;

    /** Tells whether an item is in the Bloom filter
     * \param[in] item : item to test.
     * \return the presence or not of the item
     */
	virtual bool contains (const Item& item) = 0;

    /** Tells whether the 4 neighbors of the given item are in the Bloom filter.
     * The 4 neighbors are computed from the given item by adding successively
     * nucleotides 'A', 'C', 'T' and 'G'
     * Note: some implementation may not provide this service.
     * \param[in] item : starting item from which neighbors are computed.
     * \param[in] right : if true, successors are computed, otherwise predecessors are computed.
     * \return a bitset with a boolean for the 'contains' status of each neighbor
     */
	virtual std::bitset<4> contains4 (const Item& item, bool right) = 0;

    /** Tells whether the 8 neighbors of the given item are in the Bloom filter.
     * A default implementation may be two calls to contains4 with right==true and
     * right==left and by concatenating the two bitsets.
     * Note: some implementation may not provide this service.
     * \param[in] item : starting item from which neighbors are computed.
     * \return a bitset with a boolean for the 'contains' status of each neighbor
     */
    virtual std::bitset<8> contains8 (const Item& item) = 0;

    /** Get the name of the implementation class.
     * \return the class name. */
    virtual std::string  getName   () const  = 0;

    /** Return the number of 1's in the Bloom (nibble by nibble)
     * \return the weight of the Bloom filter */
    virtual unsigned long  weight () = 0;
};

/********************************************************************************/

/** \brief Bloom filter implementation
 *
 * This implementation is an abstract factorization for subclasses. It factorizes
 * a few methods and attributes.
 */
template <typename Item> class BloomContainer : public IBloom<Item>
{
public:

    /** Constructor.
     * \param[in] tai_bloom : size (in bits) of the bloom filter.
     * \param[in] nbHash : number of hash functions to use */
    BloomContainer (u_int64_t tai_bloom, size_t nbHash = 4)
        : _hash(nbHash), n_hash_func(nbHash), blooma(0), tai(tai_bloom), nchar(0), isSizePowOf2(false)
    {
        nchar  = (1+tai/8LL);
        blooma = (unsigned char *) MALLOC (nchar*sizeof(unsigned char)); // 1 bit per elem
        system::impl::System::memory().memset (blooma, 0, nchar*sizeof(unsigned char));

        /** We look whether the provided size is a power of 2 or not.
         *   => if we have a power of two, we can optimize the modulo operations. */
        isSizePowOf2 = (tai && !(tai & (tai - 1)));

        /** In case we have a power of 2^N, we set the size to 2^N-1 and use the classical trick:
         *     a % 2^N  <=>  a & (2^N-1)
         * */
        if (isSizePowOf2)  {  tai --;  }
    }

    /** Destructor. */
    virtual ~BloomContainer ()
    {
        system::impl::System::memory().free (blooma);
    }

    /** \copydoc IBloom::getNbHash */
    size_t getNbHash () const { return n_hash_func; }

    /** \copydoc Container::contains. */
    bool contains (const Item& item)
    {
        if (isSizePowOf2)
        {
            for (size_t i=0; i<n_hash_func; i++)
            {
                u_int64_t h1 = _hash (item,i) & tai;
               // if ((blooma[h1 >> 3 ] & bit_mask[h1 & 7]) != bit_mask[h1 & 7])  {  return false;  }
				if ((blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0)  {  return false;  }

            }
        }
        else
        {
            for (size_t i=0; i<n_hash_func; i++)
            {
                u_int64_t h1 = _hash (item,i) % tai;
               // if ((blooma[h1 >> 3 ] & bit_mask[h1 & 7]) != bit_mask[h1 & 7])  {  return false;  }
				if ((blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0)  {  return false;  }

            }
        }
        return true;
    }

    /** \copydoc IBloom::contains4. */
	virtual std::bitset<4> contains4 (const Item& item, bool right)
    {   throw system::ExceptionNotImplemented ();  }

    /** \copydoc IBloom::contains8. */
    virtual std::bitset<8> contains8 (const Item& item)
    {   throw system::ExceptionNotImplemented ();  }

    /** \copydoc IBloom::getArray. */
    virtual u_int8_t*& getArray     ()  { return blooma; }

    /** \copydoc IBloom::getSize. */
    virtual u_int64_t  getSize      ()  { return nchar;  }

    /** \copydoc IBloom::getBitSize. */
    virtual u_int64_t  getBitSize   ()  { return tai;    }

    /** \copydoc IBloom::getName. */
    virtual std::string  getName    () const  = 0;

protected:

    HashFunctors<Item> _hash;
    size_t n_hash_func;

    u_int8_t* blooma;
    u_int64_t tai;
    u_int64_t nchar;
    bool      isSizePowOf2;
};

/********************************************************************************/
/** \brief Bloom filter implementation
 */
template <typename Item> class Bloom : public BloomContainer<Item>
{
public:

    /** \copydoc BloomContainer::BloomContainer */
    Bloom (u_int64_t tai_bloom, size_t nbHash = 4)  : BloomContainer<Item> (tai_bloom, nbHash)  {}

    /** \copydoc Bag::insert. */
    void insert (const Item& item)
    {
        if (this->isSizePowOf2)
        {
            for (size_t i=0; i<this->n_hash_func; i++)
            {
                u_int64_t h1 = this->_hash (item,i) & this->tai;
                this->blooma [h1 >> 3] |= bit_mask[h1 & 7];
            }
        }
        else
        {
            for (size_t i=0; i<this->n_hash_func; i++)
            {
                u_int64_t h1 = this->_hash (item,i) % this->tai;
                this->blooma [h1 >> 3] |= bit_mask[h1 & 7];
            }
        }
    }

    /** \copydoc Bag::flush */
    void flush ()  {}

    /** \copydoc IBloom::getName */
    std::string  getName () const { return "Bloom"; }

    /** Dump the Bloom filter bitset into a file.
     * OBSOLETE
     * \param[in] filename : file where to dump the bitset. */
    void dump (const char* filename)
    {
        FILE* file = fopen(filename,"wb");
        fwrite (this->blooma, sizeof(unsigned char), this->nchar, file);
        fclose (file);
    }

    /** \copydoc IBloom::weight */
    unsigned long weight()
    {
        const unsigned char oneBits[] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};
        long weight = 0;
        for(uint64_t index = 0; index < this->nchar; index++)
        {
            unsigned char current_char = this->blooma[index];
            weight += oneBits[current_char&0x0f];
            weight += oneBits[current_char>>4];
        }
        return weight;
    }
};

/********************************************************************************/
/** \brief Bloom filter null implementation
 *
 * This is a "null" implementation of the IBloom interface. It means that all
 * the contains-like requests will return false.
 */
template <typename Item> class BloomNull : public IBloom<Item>
{
public:

    /** Destructor. */
    virtual ~BloomNull() {}

    /** \copydoc IBloom::getArray */
    u_int8_t*& getArray    () { return a; }

    /** \copydoc IBloom::getSize */
    u_int64_t  getSize     () { return 0; }

    /** \copydoc IBloom::getBitSize */
    u_int64_t  getBitSize  () { return 0; }

    /** \copydoc IBloom::getNbHash */
    size_t     getNbHash   () const { return 0; }

    /** \copydoc IBloom::getName */
    virtual std::string  getName   () const  { return "BloomNull"; }

    /** \copydoc IBloom::contains4 */
    std::bitset<4> contains4 (const Item& item, bool right)  {return std::bitset<4>();}

    /** \copydoc IBloom::contains8 */
	std::bitset<8> contains8 (const Item& item)  { return std::bitset<8>(); }

    /** \copydoc IBloom::contains */
    bool contains (const Item& item) { return false; }

    /** \copydoc IBloom::insert */
    void insert (const Item& item) {}

    /** \copydoc IBloom::flush  */
    void flush ()  {}

    /** \copydoc IBloom::weight*/
    unsigned long weight ()  { return 0;}

private:
    u_int8_t* a;
};

/********************************************************************************/
/** \brief Bloom filter implementation with synchronization
 *
 * This implementation allows to insert items in the same Bloom filter by
 * different threads at the same time. It uses low level synchronization
 * mechanism like __sync_fetch_and_or
 */
template <typename Item> class BloomSynchronized : public Bloom<Item>
{
public:

    /** \copydoc Bloom::Bloom */
    BloomSynchronized (u_int64_t tai_bloom, size_t nbHash = 4)  : Bloom<Item> (tai_bloom, nbHash)  {}

    /** \copydoc Bag::insert. */
    void insert (const Item& item)
    {
        if (this->isSizePowOf2)
        {
            for (size_t i=0; i<this->n_hash_func; i++)
            {
                u_int64_t h1 = this->_hash (item,i) & this->tai;
                __sync_fetch_and_or (this->blooma + (h1 >> 3), bit_mask[h1 & 7]);
            }
        }
        else
        {
            for (size_t i=0; i<this->n_hash_func; i++)
            {
                u_int64_t h1 = this->_hash (item,i) % this->tai;
                __sync_fetch_and_or (this->blooma + (h1 >> 3), bit_mask[h1 & 7]);
            }
        }
    }

    /** \copydoc IBloom::getName*/
    std::string  getName () const { return "basic"; }
};

/********************************************************************************/
/** \brief Bloom filter implementation with CPU cache consideration
 *
 * This implementation tries to avoid CPU cache misses by improving memory locality.
 *
 * The idea is to compute the first hash function as usual. Then the other hash
 * functions are computed in order to return values near to the first value.
 *
 * The proximity is defined by a block size. Note that a too small block will
 * produce more false positive than usual.
 */
template <typename Item> class BloomCacheCoherent : public Bloom<Item>
{
public:
    
    /** Constructor.
     * \param[in] tai_bloom : size (in bits) of the bloom filter.
     * \param[in] nbHash : number of hash functions to use
     * \param[in] block_nbits : size of the block (actual 2^nbits) */
    BloomCacheCoherent (u_int64_t tai_bloom, size_t nbHash = 4,size_t block_nbits = 12)
        : Bloom<Item> (tai_bloom + 2*(1<<block_nbits), nbHash),_nbits_BlockSize(block_nbits)
    {
        _mask_block = (1<<_nbits_BlockSize) - 1;
        _reduced_tai = this->tai -  2*(1<<_nbits_BlockSize) ;//2* for neighbor coherent
    }
    
     /** \copydoc Bag::insert. */
    void insert (const Item& item)
    {
        //for insert, no prefetch, perf is not important
        u_int64_t h0;

        h0 = this->_hash (item,0) % _reduced_tai;

        __sync_fetch_and_or (this->blooma + (h0 >> 3), bit_mask[h0 & 7]);

        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 = h0  + (simplehash16( item, i) & _mask_block )   ;
            __sync_fetch_and_or (this->blooma + (h1 >> 3), bit_mask[h1 & 7]);
        }
    }
    
    /** \copydoc IBloom::getName*/
    std::string  getName () const { return "cache"; }

    /** \copydoc IBloom::getBitSize*/
    u_int64_t  getBitSize   ()  { return _reduced_tai;    }
        
    /** \copydoc Container::contains. */
    bool contains (const Item& item)
    {
        u_int64_t tab_keys [20];
        u_int64_t h0;

        h0 = this->_hash (item,0) % _reduced_tai;
        __builtin_prefetch(&(this->blooma [h0 >> 3] ), 0, 3); //preparing for read

        //compute all hashes during prefetch
        for (size_t i=1; i<this->n_hash_func; i++)
        {
           tab_keys[i] =  h0  + (simplehash16( item, i) & _mask_block );// with simplest hash
        }

        if ((this->blooma[h0 >> 3 ] & bit_mask[h0 & 7]) ==0 )  {  return false;  }
        
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 = tab_keys[i];
			if ((this->blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0)  {  return false;  }
        }
        return true;
    }
    
    /** \copydoc IBloom::weight*/
    unsigned long weight()
    {
        throw system::ExceptionNotImplemented();
    }
    
protected:
    u_int64_t _mask_block;
    size_t    _nbits_BlockSize;
    u_int64_t _reduced_tai;
};
	
/********************************************************************************/

/** \brief Bloom filter implementation with cache consideration
 *
 * This implementation improve memory locality in the Bloom filter between a kmer
 * and its neighbors. This means that this implementation should be used only
 * with Item being a kmer.
 *
 * In particular, it implements contains4 and contains8 in a clever way.
 */
template <typename Item> class BloomNeighborCoherent : public BloomCacheCoherent<Item>
{
public:

    /** Constructor.
     * \param[in] tai_bloom : size (in bits) of the bloom filter.
     * \param[in] kmersize : kmer size
     * \param[in] nbHash : number of hash functions to use
     * \param[in] block_nbits : size of the block (actual 2^nbits) */
    BloomNeighborCoherent (u_int64_t tai_bloom, size_t kmersize , size_t nbHash = 4,size_t block_nbits = 12 )  :
    BloomCacheCoherent<Item> (tai_bloom , nbHash,block_nbits), _kmerSize(kmersize)
    {
        cano2[ 0] = 0;
        cano2[ 1] = 1;
        cano2[ 2] = 2;
        cano2[ 3] = 3;
        cano2[ 4] = 4;
        cano2[ 5] = 5;
        cano2[ 6] = 3;
        cano2[ 7] = 7;
        cano2[ 8] = 8;
        cano2[ 9] = 9;
        cano2[10] = 0;
        cano2[11] = 4;
        cano2[12] = 9;
        cano2[13] = 13;
        cano2[14] = 1;
        cano2[15] = 5;

        Item un;
        un.setVal(1);
        _maskkm2  = (un << ((_kmerSize-2)*2)) - un;
        _kmerMask = (un << (_kmerSize*2))     - un;

        Item trois;
        trois.setVal(3);

        _prefmask = trois << ((_kmerSize-1)*2); //bug was here 3 instead of item trois
    }

    /** \copydoc Bag::insert. */
    void insert (const Item& item)
    {
        u_int64_t h0;
        u_int64_t racine;

        Item suffix = item & 3 ;
        Item prefix = (item & _prefmask)  >> ((_kmerSize-2)*2);
        prefix += suffix;
        prefix = prefix  & 15 ;

        u_int64_t pref_val = cano2[prefix.getVal()]; //get canonical of pref+suffix

        Item hashpart = ( item >> 2 ) & _maskkm2 ;  // delete 1 nt at each side
        Item rev =  revcomp(hashpart,_kmerSize-2);
        if(rev<hashpart) hashpart = rev; //transform to canonical

        // Item km = item;
        // rev =  revcomp(km,_kmerSize);
        // if(rev < km) km = rev; //transform to canonical
        
        racine = ((this->_hash (hashpart,0) ) % this->_reduced_tai) ;
        //h0 = ((this->_hash (item >> 2,0) ) % this->_reduced_tai)  + (suffix_val & this->_mask_block);
        //h0 = racine + (this->_hash (km,0)  & this->_mask_block);
        h0 = racine + (pref_val );
        __sync_fetch_and_or (this->blooma + (h0 >> 3), bit_mask[h0 & 7]);

        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 = h0  + ( (simplehash16( hashpart, i))  & this->_mask_block )   ;
            //	u_int64_t h1 = racine  + ( (simplehash16( km, i))  & this->_mask_block )   ; //ceci avec simplehash+8  semble ok
            //	u_int64_t h1 = h0  +  ( (this->_hash (item>>2,i)+ suffix_val)  & _mask_block );
            __sync_fetch_and_or (this->blooma + (h1 >> 3), bit_mask[h1 & 7]);
        }
    }

    /** \copydoc IBloom::getName*/
    std::string  getName () const { return "neighbor"; }

    /** \copydoc IBloom::getBitSize*/
    u_int64_t  getBitSize   ()  { return this->_reduced_tai;    }

    /** \copydoc Container::contains. */
    bool contains (const Item& item)
    {
        u_int64_t racine;

        Item suffix = item & 3 ;
        Item prefix = (item & _prefmask)  >> ((_kmerSize-2)*2);
        prefix += suffix;
        prefix = prefix  & 15 ;

        u_int64_t pref_val = cano2[prefix.getVal()]; //get canonical of pref+suffix

        Item hashpart = ( item >> 2 ) & _maskkm2 ;  // delete 1 nt at each side
        Item rev =  revcomp(hashpart,_kmerSize-2);
        if(rev<hashpart) hashpart = rev; //transform to canonical

        // Item km = item;
        // rev =  revcomp(km,_kmerSize);
        // if(rev < km) km = rev; //transform to canonical

        u_int64_t tab_keys [20];
        u_int64_t h0;

        racine = ((this->_hash (hashpart,0) ) % this->_reduced_tai) ;
        //h0 = racine + (this->_hash (km,0)  & this->_mask_block);
        h0 = racine + (pref_val  );
        //printf("h0 %llu\n",h0);

        __builtin_prefetch(&(this->blooma [h0 >> 3] ), 0, 3); //preparing for read

        //compute all hashes during prefetch
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            tab_keys[i] =  h0  + (  (simplehash16( hashpart, i)  ) & this->_mask_block );// with simplest hash
            // tab_keys[i] =  racine  + (  (simplehash16( km, i)  ) & this->_mask_block );
            // tab_keys[i] =  h0  + (  (this->_hash (item>>2,i) + suffix_val ) & _mask_block );// with simplest hash
        }

        if ((this->blooma[h0 >> 3 ] & bit_mask[h0 & 7]) == 0 )  {  return false;  } //was != bit_mask[h0 & 7]

        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 = tab_keys[i];
            if ((this->blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0 )  {  return false;  } //was != bit_mask[h0 & 7]
        }
        return true;
    }

    /** \copydoc IBloom::contains4*/
    std::bitset<4> contains4 (const Item& item, bool right)
    {
        //ask for all 4 neighbors of item in one call
        // give it CAAA
        // if right == true  it wil test the 4 neighbors AAAA, AAAC, AAAT, AAAG
        // if right == false : left extension   ACAA, CCAA , TCAA  , GCAA

        u_int64_t h0, i0, j0, k0;
        Item elem,hashpart,rev ;
        Item un ; un.setVal(1);
        Item deux ; deux.setVal(2);
        Item trois ; trois.setVal(3);

        size_t shifts = (_kmerSize -1)*2;

        if (right)  {  elem = (item << 2) & _kmerMask ;  }
        else        {  elem = (item >> 2) ;              }

        //get the canonical of middle part
        hashpart = ( elem >> 2 ) & _maskkm2 ;
        rev =  revcomp(hashpart,_kmerSize-2);
        if(rev<hashpart) hashpart = rev;

        u_int64_t racine = ((this->_hash (hashpart,0) ) % this->_reduced_tai) ;

        __builtin_prefetch(&(this->blooma [racine >> 3] ), 0, 3); //preparing for read

        Item tmp,suffix,prefix;
        u_int64_t pref_val;

        //with val of prefix+suffix  for neighbor shift

        tmp = elem;
        suffix = tmp & 3 ;
        prefix = (tmp & _prefmask)  >> ((_kmerSize-2)*2);
        prefix += suffix;
        pref_val = cano2[prefix.getVal()]; //get canonical of pref+suffix

        h0 = racine + (pref_val  & this->_mask_block);

        if(right) tmp = elem+un;
        else tmp = elem + (un<<shifts) ;
        suffix = tmp & 3 ;
        prefix = (tmp & _prefmask)  >> ((_kmerSize-2)*2);
        prefix += suffix;
        pref_val = cano2[prefix.getVal()]; //get canonical of pref+suffix

        i0 = racine + (pref_val  & this->_mask_block);

        if(right) tmp = elem+deux;
        else tmp = elem + (deux<<shifts) ;
        suffix = tmp & 3 ;
        prefix = (tmp & _prefmask)  >> ((_kmerSize-2)*2);
        prefix += suffix;
        pref_val = cano2[prefix.getVal()]; //get canonical of pref+suffix

        j0 = racine + (pref_val  & this->_mask_block);

        if(right) tmp = elem+trois;
        else tmp = elem + (trois<<shifts) ;
        suffix = tmp & 3 ;
        prefix = (tmp & _prefmask)  >> ((_kmerSize-2)*2);
        prefix += suffix;
        pref_val = cano2[prefix.getVal()]; //get canonical of pref+suffix

        k0 = racine + (pref_val  & this->_mask_block);

        /*
         //with full hash of kmer for neighbor shift
        tmp = elem;
        rev =  revcomp(tmp,_kmerSize);
        if(rev < tmp) tmp = rev;
        h0 = racine + (this->_hash (tmp,0)  & this->_mask_block);

        if(right) tmp = elem+un;
        else tmp = elem + (un<<shifts) ;
        rev =  revcomp(tmp,_kmerSize); // many revcomp, optim possible
        if(rev < tmp) tmp = rev;
        i0 = racine + (this->_hash (tmp,0)  & this->_mask_block);

        if(right) tmp = elem+deux;
        else tmp = elem + (deux<<shifts) ;
        rev =  revcomp(tmp,_kmerSize);
        if(rev < tmp) tmp = rev;
        j0 = racine + (this->_hash (tmp,0)  & this->_mask_block);

        if(right) tmp = elem+trois;
        else tmp = elem + (trois<<shifts) ;
        rev =  revcomp(tmp,_kmerSize);
        if(rev < tmp) tmp = rev;
        k0 = racine + (this->_hash (tmp,0)  & this->_mask_block);
         */

        u_int64_t tab_hashes [20];

        //compute all hashes during prefetch
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            tab_hashes[i] = simplehash16( hashpart, i) & this->_mask_block ;
        }

        std::bitset<4> resu;
        resu.set (0, true);
        resu.set (1, true);
        resu.set (2, true);
        resu.set (3, true);

        if ((this->blooma[h0 >> 3 ] & bit_mask[h0 & 7]) == 0)  {  resu.set (0, false);  }
        if ((this->blooma[i0 >> 3 ] & bit_mask[i0 & 7]) == 0)  {  resu.set (1, false);  }
        if ((this->blooma[j0 >> 3 ] & bit_mask[j0 & 7]) == 0)  {  resu.set (2, false);  }
        if ((this->blooma[k0 >> 3 ] & bit_mask[k0 & 7]) == 0)  {  resu.set (3, false);  }

        //plus rapide avec 4 boucles separees que une ci dessous avec test pour break
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 =  h0 +  tab_hashes[i]   ;
            if ( (this->blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0  )  {  resu.set (0, false); break; }
        }
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 =  i0 +  tab_hashes[i]   ;
            if ( (this->blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0  )  {  resu.set (1, false); break; }
        }
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 =  j0 +  tab_hashes[i]   ;
            if ( (this->blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0  )  {  resu.set (2, false); break; }
        }
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 =  k0 +  tab_hashes[i]   ;
            if ( (this->blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0  )  {  resu.set (3, false); break; }
        }

        /*
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 =  h0 +  tab_hashes[i]   ;
            u_int64_t i1 =  i0 +  tab_hashes[i]   ;
            u_int64_t j1 =  j0 +  tab_hashes[i]   ;
            u_int64_t k1 =  k0 +  tab_hashes[i]   ;

            if (resu[0] && (this->blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0  )  {  resu[0]=false;  }  //test  resu[0] &&
            if (resu[1] &&(this->blooma[i1 >> 3 ] & bit_mask[i1 & 7]) == 0  )  {  resu[1]=false;  }
            if (resu[2] &&(this->blooma[j1 >> 3 ] & bit_mask[j1 & 7]) == 0  )  {  resu[2]=false;  }
            if (resu[3] &&(this->blooma[k1 >> 3 ] & bit_mask[k1 & 7]) == 0  )  {  resu[3]=false;  }

            if(resu[0]== false &&  resu[1]== false && resu[2]== false && resu[3]== false)
                break;
        }
         */

        return resu;
    }

    /** \copydoc IBloom::contains8*/
    std::bitset<8> contains8 (const Item& item)
    {
        std::bitset<4> resultRight = this->contains4 (item, true);
        std::bitset<4> resultLeft  = this->contains4 (item, false);
        std::bitset<8> result;
        size_t i=0;
        for (size_t j=0; j<4; j++)  { result.set (i++, resultRight[j]); }
        for (size_t j=0; j<4; j++)  { result.set (i++, resultLeft [j]); }
        return result;
    }

private:
    unsigned int cano2[16];
    Item _maskkm2;
    Item _prefmask;
    Item _kmerMask;
    size_t _kmerSize;
};
    
/********************************************************************************/

/** \brief Bloom filter implementation with cache consideration
 *
 * This implementation improve memory locality in the Bloom filter between a kmer
 * and its neighbors. This means that this implementation should be used only
 * with Item being a kmer.
 *
 * In particular, it implements contains4 and contains8 in a clever way.
 */
template <typename Item> class BloomExtendedNeighborCoherent : public BloomCacheCoherent<Item>
{
public:

    /** Constructor.
     * \param[in] tai_bloom : size (in bits) of the bloom filter.
     * \param[in] kmersize : kmer size
     * \param[in] hmersize : hashpart size
     * \param[in] nbHash : number of hash functions to use
     * \param[in] block_nbits : size of the block (actual 2^nbits) */
    BloomExtendedNeighborCoherent (u_int64_t tai_bloom, size_t kmersize, size_t nbHash = 7, size_t block_nbits = 12)
        : BloomCacheCoherent<Item> (tai_bloom, nbHash,block_nbits), _kmerSize(kmersize)
    {
        _smerSize = _kmerSize - 2;
        _hmerSize = _smerSize - 8;

        cano6 = (unsigned short int *) MALLOC (0x1000 * sizeof(unsigned short int));
        system::impl::System::memory().memset (cano6, 0, 0x1000 * sizeof(unsigned short int));

        hpos = (unsigned char *) MALLOC (0x40000 * sizeof(unsigned char));
        system::impl::System::memory().memset (hpos, 0, 0x40000 * sizeof(unsigned char));

        precomputeCano6();
        precomputeHpos();

        Item un = 1;
        _kmerMask = (un << (_kmerSize*2))   - un;
        _smerMask = (un << (_smerSize*2))   - un;
        _hmerMask = (un << (_hmerSize*2))   - un;

        Item trois = 3;
        _kmerPrefMask = ((Item)0x3f) << ((_kmerSize-3)*2);
        _smerPrefMask = trois << ((_smerSize-1)*2);

        _hmerCount = _smerSize - _hmerSize + 1;

        _sharedpart = _smerMask + un; // > max value
        _hashpartFwd = _hmerMask + un; // > max value
        _hashpartRev = _hmerMask + un; // > max value

        _hashpartHits = 0;
    }

    ~BloomExtendedNeighborCoherent() {
        system::impl::System::memory().free (cano6);
        system::impl::System::memory().free (hpos);
    }

    /** \copydoc Bag::insert. */
    void insert (const Item& item)
    {
        u_int64_t h0;
        u_int64_t racine;

        Item suffix = item & ((Item)0x3f);
        Item limits = (item & _kmerPrefMask)  >> ((_kmerSize-6)*2);
        limits += suffix;
        u_int64_t delta = cano6[limits.getVal()];

        Item sharedpart = (item >> 2) & _smerMask;  // delete 1 nt at each side
        Item rev =  revcomp(sharedpart, _smerSize);
        if(rev < sharedpart) sharedpart = rev; //transform to canonical

        _hpart = extractHashpart(sharedpart);
        _hpartHash = this->_hash (_hpart, 0);
        
        racine = _hpartHash % this->_reduced_tai;
        h0 = racine + delta;
        __sync_fetch_and_or(this->blooma + (h0 >> 3), bit_mask[h0 & 7]);


        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 = h0 + ((simplehash16( _hpart, i)) & this->_mask_block);
            __sync_fetch_and_or(this->blooma + (h1 >> 3), bit_mask[h1 & 7]);
        }
    }

    /** \copydoc IBloom::getName*/
    std::string  getName () const { return "neighbor2"; }

    /** \copydoc IBloom::getBitSize*/
    u_int64_t  getBitSize   ()  { return this->_reduced_tai;    }

    /** \copydoc Container::contains. */
    bool contains (const Item& item, const Item& next = 0)
    {
        Item suffix = item & ((Item)0x3f);
        Item limits = (item & _kmerPrefMask)  >> ((_kmerSize-6)*2);
        limits += suffix;
        u_int64_t delta = cano6[limits.getVal()];

        bool reverse = false;
        Item sharedpart = (item >> 2) & _smerMask ;  // delete 1 nt at each side
        Item rev =  revcomp(sharedpart, _smerSize);
        if(rev < sharedpart) {
            sharedpart = rev; //transform to canonical
            reverse = true;
        }

        if (_sharedpart != sharedpart)
        {
            Item hpart = extractHashpart(sharedpart);

            if (reverse)
            {
                if (hpart != _hashpartRev) {
                    _hashpartRev = hpart;
                    _hashpartRevHash = this->_hash(_hashpartRev, 0);
                    prepareTabHashes(_tabHashesRev, _hashpartRev);
                }
                else
                {
                    _hashpartHits++;
                }

                _hpart = reverse ? _hashpartRev : _hashpartFwd;
                _hpartHash = _hashpartRevHash;
                _tabHashes = _tabHashesRev;
            }
            else
            {
                 if (hpart != _hashpartFwd) {
                    _hashpartFwd = hpart;
                    _hashpartFwdHash = this->_hash(_hashpartFwd, 0);
                    prepareTabHashes(_tabHashesFwd, _hashpartFwd);
                }
                else
                {
                    _hashpartHits++;
                }

                _hpart = _hashpartFwd;
                _hpartHash = _hashpartFwdHash;
                _tabHashes = _tabHashesFwd;
            }

            _sharedpart = sharedpart;
        }
        
        u_int64_t racine, h0;

        racine = _hpartHash % this->_reduced_tai;
        h0 = racine + delta;

        __builtin_prefetch(&(this->blooma [racine >> 3] ), 0, 3); // preparing for read

        // compute all hashes during prefetch
        u_int64_t tab_keys [20];
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            tab_keys[i] =  h0 + _tabHashes[i]; // with simplest hash
        }

        if ((this->blooma[h0 >> 3 ] & bit_mask[h0 & 7]) == 0 )  {  return false;  } // was != bit_mask[h0 & 7]

        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 = tab_keys[i];
            if ((this->blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0 )  {  return false;  } // was != bit_mask[h1 & 7]
        }

        return true;
    }

    /** \copydoc IBloom::contains4*/
    std::bitset<4> contains4 (const Item& item, bool right)
    {
        //ask for all 4 neighbors of item in one call
        // give it CAAA
        // if right == true  it wil test the 4 neighbors AAAA, AAAC, AAAT, AAAG
        // if right == false : left extension   ACAA, CCAA , TCAA  , GCAA

        u_int64_t h0, i0, j0, k0;
        Item elem,sharedpart,rev ;
        Item un = 1;
        Item deux = 2;
        Item trois = 3;
        // soleil !
        bool reverse = false;

        size_t shifts = (_kmerSize -1)*2;

        if (right)  {  elem = (item << 2) & _kmerMask ;  }
        else        {  elem = (item >> 2) ;              }

        //get the canonical of middle part
        sharedpart = ( elem >> 2 ) & _smerMask ;
        rev =  revcomp(sharedpart,_smerSize);
        if(rev<sharedpart) {
            sharedpart = rev;
            reverse = true;
        }


            Item hpart = extractHashpart(sharedpart);

            if (reverse)
            {
                if (hpart != _hashpartRev) {
                    _hashpartRev = hpart;
                    _hashpartRevHash = this->_hash(_hashpartRev, 0);
                    prepareTabHashes(_tabHashesRev, _hashpartRev);
                }
                else
                {
                    _hashpartHits++;
                }

                _hpart = reverse ? _hashpartRev : _hashpartFwd;
                _hpartHash = _hashpartRevHash;
                _tabHashes = _tabHashesRev;
            }
            else
            {
                 if (hpart != _hashpartFwd) {
                    _hashpartFwd = hpart;
                    _hashpartFwdHash = this->_hash(_hashpartFwd, 0);
                    prepareTabHashes(_tabHashesFwd, _hashpartFwd);
                }
                else
                {
                    _hashpartHits++;
                }

                _hpart = _hashpartFwd;
                _hpartHash = _hashpartFwdHash;
                _tabHashes = _tabHashesFwd;
            }

            _sharedpart = sharedpart;


        u_int64_t racine = _hpartHash % this->_reduced_tai;

        __builtin_prefetch(&(this->blooma [racine >> 3] ), 0, 3); //preparing for read

        Item tmp,suffix,limits;
        u_int64_t delta;

        //with val of prefix+suffix  for neighbor shift

        tmp = elem;
        suffix = tmp & ((Item)0x3f) ;
        limits = (tmp & _kmerPrefMask)  >> ((_kmerSize-6)*2);
        limits += suffix;
        delta = cano6[limits.getVal()]; //get canonical of pref+suffix

        h0 = racine + delta;

        if(right) tmp = elem+un;
        else tmp = elem + (un<<shifts) ;
        suffix = tmp & ((Item)0x3f) ;
        limits = (tmp & _kmerPrefMask)  >> ((_kmerSize-6)*2);
        limits += suffix;
        delta = cano6[limits.getVal()]; //get canonical of pref+suffix

        i0 = racine + delta;

        if(right) tmp = elem+deux;
        else tmp = elem + (deux<<shifts) ;
        suffix = tmp & ((Item)0x3f) ;
        limits = (tmp & _kmerPrefMask)  >> ((_kmerSize-6)*2);
        limits += suffix;
        delta = cano6[limits.getVal()]; //get canonical of pref+suffix

        j0 = racine + delta;

        if(right) tmp = elem+trois;
        else tmp = elem + (trois<<shifts) ;
        suffix = tmp & ((Item)0x3f) ;
        limits = (tmp & _kmerPrefMask)  >> ((_kmerSize-6)*2);
        limits += suffix;
        delta = cano6[limits.getVal()]; //get canonical of pref+suffix

        k0 = racine + delta;


        std::bitset<4> resu;
        resu.set (0, true);
        resu.set (1, true);
        resu.set (2, true);
        resu.set (3, true);

        if ((this->blooma[h0 >> 3 ] & bit_mask[h0 & 7]) == 0)  {  resu.set (0, false);  }
        if ((this->blooma[i0 >> 3 ] & bit_mask[i0 & 7]) == 0)  {  resu.set (1, false);  }
        if ((this->blooma[j0 >> 3 ] & bit_mask[j0 & 7]) == 0)  {  resu.set (2, false);  }
        if ((this->blooma[k0 >> 3 ] & bit_mask[k0 & 7]) == 0)  {  resu.set (3, false);  }

        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t h1 =  h0 +  _tabHashes[i]   ;
            if ( (this->blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0  )  {  resu.set (0, false); break; }
        }
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t i1 =  i0 +  _tabHashes[i]   ;
            if ( (this->blooma[i1 >> 3 ] & bit_mask[i1 & 7]) == 0  )  {  resu.set (1, false); break; }
        }
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t j1 =  j0 +  _tabHashes[i]   ;
            if ( (this->blooma[j1 >> 3 ] & bit_mask[j1 & 7]) == 0  )  {  resu.set (2, false); break; }
        }
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            u_int64_t k1 =  k0 +  _tabHashes[i]   ;
            if ( (this->blooma[k1 >> 3 ] & bit_mask[k1 & 7]) == 0  )  {  resu.set (3, false); break; }
        }

        return resu;
    }

    u_int64_t getHashpartHits() const {
        return _hashpartHits;
    }


private:
    unsigned short int *cano6;
    unsigned char *hpos;

    Item _kmerMask;
    Item _smerMask;
    Item _hmerMask;
    Item _kmerPrefMask;
    Item _smerPrefMask;
    size_t _kmerSize;
    size_t _smerSize;
    size_t _hmerSize;
    int16_t _hmerCount;

    Item _sharedpart;

    Item _hpart;
    u_int64_t _hpartHash;
    u_int64_t* _tabHashes;

    Item _hashpartFwd;
    u_int64_t _hashpartFwdHash;
    u_int64_t _tabHashesFwd[20];

    Item _hashpartRev;
    u_int64_t _hashpartRevHash;
    bool _hpartRevHashComputed;
    u_int64_t _tabHashesRev[20];

    u_int64_t _hashpartHits;


    inline void prepareTabHashes(u_int64_t* tab, const Item& toHash)
    {
        for (size_t i=1; i<this->n_hash_func; i++)
        {
            tab[i] = simplehash16(toHash, i) & this->_mask_block ;
        }
    }

    Item extractHashpart(const Item& sharedpart)
    {
        Item posPart = sharedpart >> (_smerSize*2 - 18);
        unsigned char pos = hpos[posPart.getVal()];

        Item hpart = (sharedpart >> ((_hmerCount - pos - 1)*2)) & _hmerMask;
        return hpart;
    }

    void precomputeCano6()
    {
        for (uint64_t i=0; i<0x1000; i++) {
            Item cur = (Item)i;
            Item rev = revcomp(cur, 6);
            Item cano = (cur < rev) ? cur : rev;
            cano6[i] = cano.getVal();
        }
    }

    unsigned char minpos(const u_int64_t& nmer, size_t n)
    {
        static const size_t MINIMIZER_SIZE = 2;
        static const uint64_t minMask = (1 << (MINIMIZER_SIZE*2)) - 1;

        uint64_t min = nmer & minMask;
        unsigned char pos = n - MINIMIZER_SIZE;

        for (size_t i=1; i <= n - MINIMIZER_SIZE; i++) {
            uint64_t cur = (nmer >> (i*2)) & minMask;

            if (cur < min) {
                min = cur;
                pos = n - MINIMIZER_SIZE - i;
            }
        }

        return pos;
    }

    void precomputeHpos()
    {
        for (uint64_t i=0; i<0x40000; i++) {
            hpos[i] = minpos(i, 9);
        }
    }
};
	
/********************************************************************************/

/** \brief Factory that creates IBloom instances
 *
 */
class BloomFactory
{
public:

    /** Singleton method
     * \return the singleton. */
    static BloomFactory& singleton()  { static BloomFactory instance; return instance; }

    /** Create a IBloom instance
     * \param[in] kind : kind of the IBloom instance to be created
     * \param[in] tai_bloom : size of the Bloom filter (in bits)
     * \param[in] nbHash : number of hash functions for the Bloom filter
     * \param[in] kmersize : kmer size (used only for some implementations).
     */
    template<typename T> IBloom<T>* createBloom (tools::misc::BloomKind kind, u_int64_t tai_bloom, size_t nbHash, size_t kmersize)
    {
        switch (kind)
        {
            case tools::misc::BLOOM_NONE:      return new BloomNull<T>             ();
            case tools::misc::BLOOM_BASIC:     return new BloomSynchronized<T>     (tai_bloom, nbHash);
            case tools::misc::BLOOM_CACHE:     return new BloomCacheCoherent<T>    (tai_bloom, nbHash);
			case tools::misc::BLOOM_NEIGHBOR:  return new BloomNeighborCoherent<T> (tai_bloom, kmersize, nbHash);
            case tools::misc::BLOOM_DEFAULT:   return new BloomCacheCoherent<T>    (tai_bloom, nbHash);
            default:        throw system::Exception ("bad Bloom kind %d in createBloom", kind);
        }
    }

    /** Create a IBloom instance
     * \param[in] name : kind name of the IBloom instance to be created
     * \param[in] sizeStr : size of the Bloom filter (in bits) as a string
     * \param[in] nbHashStr : number of hash functions for the Bloom filter as a string
     * \param[in] kmerSizeStr : kmer size (used only for some implementations) as a string.
     */
    template<typename T> IBloom<T>* createBloom (
        const std::string& name,
        const std::string& sizeStr,
        const std::string& nbHashStr,
        const std::string& kmerSizeStr
    )
    {
        //std::cout << "custom createbloom, name=" << name << " size=" << sizeStr << " nbHash=" << nbHashStr << " k=" << kmerSizeStr << std::endl;
        tools::misc::BloomKind kind;  parse (name, kind);
        return createBloom<T> (kind, (u_int64_t)atol (sizeStr.c_str()), (size_t)atol (nbHashStr.c_str()), atol (kmerSizeStr.c_str()));
    }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BLOOM_HPP_ */
