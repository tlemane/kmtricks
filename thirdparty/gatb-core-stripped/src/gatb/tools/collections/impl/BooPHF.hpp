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

/** \file BooPHF.hpp
 *  \brief Minimal Perfect Hash Function from Guillaume
 */

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_BOOPHF_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_BOOPHF_HPP_

/********************************************************************************/

#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>

#include <BooPHF/BooPHF.h>

#include <random> // for mt19937_64

/********************************************************************************/
namespace gatb        {
namespace core        {
namespace tools       {
namespace collections {
namespace impl        {
/********************************************************************************/


typedef std::pair<u_int8_t const*, u_int8_t const*> byte_range_t;

/** For some specialization (see below), we need to adapt the key type to some
 * range of raw data in memory. We provide here a default adaptor that can
 * be used as default template type for the MPHF class.
 */
template<typename T>
struct AdaptatorDefault
{
    byte_range_t operator() (const T& t) const
    {
        const u_int8_t* buf = reinterpret_cast <u_int8_t const*> (&t);
        const u_int8_t* end = buf + sizeof(T);
        return byte_range_t(buf, end);
    }
};

// from emphf, https://github.com/ot/emphf/blob/master/base_hash.hpp
// Apache License 2
// itself was adapted from http://www.burtleburtle.net/bob/c/lookup8.c
inline uint64_t unaligned_load64(uint8_t const* from)
{
    uint64_t tmp;
    memcpy(reinterpret_cast<char*>(&tmp), from, 8);
          //(ot): reverse bytes in big-endian architectures
         return tmp;
    }
struct jenkins64_hasher {

	typedef uint64_t seed_t;
	typedef uint64_t hash_t;
	typedef std::tuple<hash_t, hash_t, hash_t> hash_triple_t;

	jenkins64_hasher()
	{}

	jenkins64_hasher(uint64_t seed)
		: m_seed(seed)
	{}

	template <typename Rng>
		static jenkins64_hasher generate(Rng& rng)
		{
			return jenkins64_hasher(rng());
		}

	// Adapted from http://www.burtleburtle.net/bob/c/lookup8.c
	hash_triple_t operator()(byte_range_t s) const
	{
		using std::get;
		hash_triple_t h(m_seed, m_seed, 0x9e3779b97f4a7c13ULL);

		size_t len = (size_t)(s.second - s.first);
		uint8_t const* cur = s.first;
		uint8_t const* end = s.second;

		while (end - cur >= 24) {
			get<0>(h) += unaligned_load64(cur);
			cur += 8;
			get<1>(h) += unaligned_load64(cur);
			cur += 8;
			get<2>(h) += unaligned_load64(cur);
			cur += 8;

			mix(h);
		}

		get<2>(h) += len;

		switch (end - cur) {
			case 23: get<2>(h) += (uint64_t(cur[22]) << 56);
			case 22: get<2>(h) += (uint64_t(cur[21]) << 48);
			case 21: get<2>(h) += (uint64_t(cur[20]) << 40);
			case 20: get<2>(h) += (uint64_t(cur[19]) << 32);
			case 19: get<2>(h) += (uint64_t(cur[18]) << 24);
			case 18: get<2>(h) += (uint64_t(cur[17]) << 16);
			case 17: get<2>(h) += (uint64_t(cur[16]) << 8);
					 // the first byte of c is reserved for the length
			case 16: get<1>(h) += (uint64_t(cur[15]) << 56);
			case 15: get<1>(h) += (uint64_t(cur[14]) << 48);
			case 14: get<1>(h) += (uint64_t(cur[13]) << 40);
			case 13: get<1>(h) += (uint64_t(cur[12]) << 32);
			case 12: get<1>(h) += (uint64_t(cur[11]) << 24);
			case 11: get<1>(h) += (uint64_t(cur[10]) << 16);
			case 10: get<1>(h) += (uint64_t(cur[ 9]) << 8);
			case  9: get<1>(h) += (uint64_t(cur[ 8]));
			case  8: get<0>(h) += (uint64_t(cur[ 7]) << 56);
			case  7: get<0>(h) += (uint64_t(cur[ 6]) << 48);
			case  6: get<0>(h) += (uint64_t(cur[ 5]) << 40);
			case  5: get<0>(h) += (uint64_t(cur[ 4]) << 32);
			case  4: get<0>(h) += (uint64_t(cur[ 3]) << 24);
			case  3: get<0>(h) += (uint64_t(cur[ 2]) << 16);
			case  2: get<0>(h) += (uint64_t(cur[ 1]) << 8);
			case  1: get<0>(h) += (uint64_t(cur[ 0]));
			case 0: break; // nothing to add
			default: assert(false);
		}

		mix(h);

		return h;
	}

	// rehash a hash triple
	hash_triple_t operator()(hash_triple_t h) const
	{
		std::get<0>(h) += m_seed;
		std::get<1>(h) += m_seed;
		std::get<2>(h) += 0x9e3779b97f4a7c13ULL;

		mix(h);

		return h;
	}

	void swap(jenkins64_hasher& other)
	{
		std::swap(m_seed, other.m_seed);
	}

	void save(std::ostream& os) const
	{
		os.write(reinterpret_cast<char const*>(&m_seed), sizeof(m_seed));
	}

	void load(std::istream& is)
	{
		is.read(reinterpret_cast<char*>(&m_seed), sizeof(m_seed));
	}

	seed_t seed() const
	{
		return m_seed;
	}

	protected:

	static void mix(hash_triple_t& h)
	{
		uint64_t& a = std::get<0>(h);
		uint64_t& b = std::get<1>(h);
		uint64_t& c = std::get<2>(h);

		a -= b; a -= c; a ^= (c >> 43);
		b -= c; b -= a; b ^= (a << 9);
		c -= a; c -= b; c ^= (b >> 8);
		a -= b; a -= c; a ^= (c >> 38);
		b -= c; b -= a; b ^= (a << 23);
		c -= a; c -= b; c ^= (b >> 5);
		a -= b; a -= c; a ^= (c >> 35);
		b -= c; b -= a; b ^= (a << 49);
		c -= a; c -= b; c ^= (b >> 11);
		a -= b; a -= c; a ^= (c >> 12);
		b -= c; b -= a; b ^= (a << 18);
		c -= a; c -= b; c ^= (b >> 22);
	}

	seed_t m_seed;
};



/** \brief Minimal Perfect Hash Function
 *
 * This is a specialization of the MPHF<Key,Adaptor,exist> class for exist=true.
 * It uses BooPHF for the implementation and is most a wrapper between BooPHF and
 * GATB-CORE concepts.
 */
/** \brief Perfect minimal hash function for a given kind of key
 *
 * This class provides an interface for getting hash codes for some key type T, which
 * can be done through the operator() method
 *
 * This class is not a classic hash feature because it hashes only a given set of T items
 * (provided as a T iterator) through its 'build' method. Once building is done, hash code
 * can be accessed through the operator()
 *
 * We propose here a default implementation that doesn't do much. The idea behind is that
 * we can specialize the class for the 'exist' template argument in order to provide a true
 * implementation (through EMPHF library for instance). If such an implementation exists,
 * the constant 'enabled' will be true, which allows to test it in the code (it is a little
 * bit better than using compilation flag).
 */

template<typename Key,typename Adaptator=AdaptatorDefault<Key>, class Progress=tools::misc::impl::ProgressNone>
class BooPHF : public system::SmartPointer
{
private:

    // a hash wrapper that calls emphf's hasher to produce, given an element, a single hash value for BooPHF
    class hasher_t
    {
        typedef jenkins64_hasher BaseHasher;
        BaseHasher emphf_hasher;
        Adaptator adaptor;
            
       public:
        hasher_t(){
            std::mt19937_64 rng(37); // deterministic seed
            emphf_hasher = BaseHasher::generate(rng);
        }

 

        uint64_t operator ()  (const Key& key, uint64_t seed = 0) const  {  
                if (seed != 0x33333333CCCCCCCCULL)
                    return std::get<0>(emphf_hasher(adaptor(key)));  
                return std::get<2>(emphf_hasher(adaptor(key)));   
                // this is a big hack, because I'm lazy. 
                // I wanted to return two different hashes depending on how boophf calls it
                // since I contrl BooPHF code's, I know it calls this function with 0x33333333CCCCCCCCULL as the second seed.
                }
    };

    typedef boomphf::mphf<  Key, hasher_t  > boophf_t;

public:

    /** Definition of a hash value. */
    typedef u_int64_t Code;

    /** Constructor. */
    BooPHF () : isBuilt(false), nbKeys(0)  {}

    /** Build the hash function from a set of items.
     * \param[in] iterable : keys iterator
     * \param[in] progress : object that listens to the event of the algorithm */
    void build (tools::collections::Iterable<Key>* iterable, int nbThreads = 1, tools::dp::IteratorListener* progress=0)
    {
        if (isBuilt==true) { throw system::Exception ("MFHP: built already done"); }

        /** We create an iterator from the iterable. */
        tools::dp::Iterator<Key>* iter = iterable->iterator();
        LOCAL (iter);

        size_t nbElts = iterable->getNbItems();

        iterator_wrapper kmers (iter); 

		bool withprogress = true;

		if (progress==0)
			withprogress = false;
		

        bphf =  boophf_t(nbElts, kmers, nbThreads, 3.0 /*much faster construction than gamma=1*/, withprogress);

        isBuilt = true;
        nbKeys  = iterable->getNbItems();
    }

    /** Returns the hash code for the given key. WARNING : default implementation here will
     * throw an exception.
     * \param[in] key : the key to be hashed
     * \return the hash value. */
    Code operator () (const Key& key)
    {
        return bphf.lookup (key);
    }

    /** Returns the number of keys.
     * \return keys number */
    size_t size() const { return bphf.nbKeys(); }

    /** Load hash function from a collection*/
    size_t load (tools::storage::impl::Group& group, const std::string& name)
    {
        /** We need an input stream for the given collection given by group/name. */
        tools::storage::impl::Storage::istream is (group, name);
		bphf =  boophf_t();
        bphf.load (is);
        return size();
    }

    /** Save hash function to a collection
     * \return the number of bytes of the saved data. */
    size_t save (tools::storage::impl::Group& group, const std::string& name)
    {
        /** We need an output stream for the given collection given by group/name. */
        tools::storage::impl::Storage::ostream os (group, name);
        bphf.save (os);
        /** We set the number of keys as an attribute of the group. */
        group.addProperty ("nb_keys", misc::impl::Stringify().format("%d",nbKeys)); // FIXME: maybe overflow here
        return os.tellp();
    }

private:

    boophf_t  bphf;
    bool      isBuilt;
    size_t    nbKeys;

private:

    class iterator_adaptator : public std::iterator<std::forward_iterator_tag, const Key>
    {
    public:
        iterator_adaptator()  : iterator(0), pos(0) {}

        iterator_adaptator(tools::dp::Iterator<Key>* iterator)  : iterator(iterator), pos(0)  {  iterator->first();  }

        Key const& operator*()  {  return iterator->item();  }

        iterator_adaptator& operator++()
        {
            iterator->next();
            pos++;
            if (iterator->isDone())
            {
                iterator = nullptr;
                pos = 0;
            }
            return *this;
        }

        friend bool operator==(iterator_adaptator const& lhs, iterator_adaptator const& rhs)
        {
            if (!lhs.iterator || !rhs.iterator)  {  if (!lhs.iterator && !rhs.iterator) {  return true; } else {  return false;  } }
            return rhs.pos == lhs.pos;
        }

        friend bool operator!=(iterator_adaptator const& lhs, iterator_adaptator const& rhs)  {  return !(lhs == rhs);  }

    private:
        tools::dp::Iterator<Key>* iterator;
        unsigned long pos;
    };

    class iterator_wrapper
    {
    public:
        iterator_wrapper (tools::dp::Iterator<Key>* iterator) : iterator(iterator) {}

        iterator_adaptator begin() const  {  return iterator_adaptator (iterator); }
        iterator_adaptator end  () const  {  return iterator_adaptator ();         }
        size_t        size () const  {  return 0;                        }

    private:
        // noncopyble // FIXME: made it copyable because boophf needed it; need to see if it's correct
        //iterator_wrapper(iterator_wrapper const&);
        //iterator_wrapper& operator=(iterator_wrapper const&);
        tools::dp::Iterator<Key>* iterator;
    };
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_BOOPHF_HPP_ */
