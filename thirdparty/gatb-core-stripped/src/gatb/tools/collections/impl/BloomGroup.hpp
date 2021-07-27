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

/** \file BloomGroup.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Bloom Group implementation
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BLOOM_GROUP_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BLOOM_GROUP_HPP_

/********************************************************************************/

#include <gatb/tools/collections/impl/Bloom.hpp>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

/* EXPERIMENTAL (not documented). */
template <typename Item, size_t prec=1> class BloomGroupOld : public system::SmartPointer
{
public:

    typedef tools::math::LargeInt<prec> Result;

    /** */
    BloomGroupOld (u_int64_t size, size_t nbHash=4)
        : _hash(nbHash), _nbHash(nbHash), _size(size), _blooma(0)
    {
        _blooma = (Result*) MALLOC (_size*sizeof(Result));
        system::impl::System::memory().memset (_blooma, 0, _size*sizeof(Result));
    }

    /** */
    BloomGroupOld (const std::string& uri)
        : _hash(0), _nbHash(0), _size(0), _blooma(0)
    {
        load (uri);
    }

    /** */
    ~BloomGroupOld ()  {  FREE (_blooma); }

    /** */
    std::string getName () const { return "BloomGroupOld"; }

    /** */
    void insert (const Item& item, size_t idx)
    {
        static const Result ONE (1);

        for (size_t i=0; i<this->_nbHash; i++)
        {
            u_int64_t h1 = this->_hash (item, i) % this->_size;
#if 1
            this->_blooma[h1] |= (ONE << idx);
#else
            this->_blooma[h1].sync_fetch_and_or (ONE << idx);
#endif
        }
    }

    /** Return the size (in bytes). */
    u_int64_t getMemSize () const { return _size*sizeof(Result); }

    /** */
    void save (const std::string& uri)
    {
        system::IFile* file = system::impl::System::file().newFile (uri, "wb+");
        if (file != 0)
        {
            /** We write the nb of hash functions. */
            file->fwrite (&_nbHash, sizeof(_nbHash), 1);

            /** We write the size of the blooms. */
            file->fwrite (&_size, sizeof(_size), 1);

            /** We write the blooms info. */
            file->fwrite (_blooma, _size*sizeof(Result), 1);

#if 1
for (size_t i=0; i<10; i++)
{
    for (size_t j=0; j<prec; j++)
    {
        printf ("%8x ", _blooma[i].value[j]);
    }
    printf ("\n");
}
#endif

            delete file;
        }
    }

    /** */
    void load (const std::string& uri)
    {
        system::IFile* file = system::impl::System::file().newFile (uri, "rb+");
        if (file != 0)
        {
            /** We read the nb of hash functions. */
            file->fread (&_nbHash, sizeof(_nbHash), 1);

            /** We read the size of the blooms. */
            file->fread (&_size, sizeof(_size), 1);

            /** We allocate the array. */
            _blooma = (Result*) MALLOC (_size*sizeof(Result));
            system::impl::System::memory().memset (_blooma, 0, _size*sizeof(Result));

            /** We read the blooms info. */
            file->fread (_blooma, _size*sizeof(Result), 1);

            delete file;
        }
    }

    /** */
    bool contains (const Item& item, size_t idx)
    {
        static const Result ONE (1);
        for (size_t i=0; i<this->_nbHash; i++)
        {
            u_int64_t h1 = this->_hash (item, i) % this->_size;
            if ( (_blooma[h1] & (ONE << idx)) != (ONE << idx) )  {  return false;  }
        }
        return true;
    }

    /** */
    Result contains (const Item& item)
    {
        static const Result ZERO (0);
        Result res = ~ZERO;

        for (size_t i=0; i<this->_nbHash; i++)
        {
            u_int64_t h1 = this->_hash (item, i) % this->_size;
            res &= _blooma [h1];
        }
        return res;
    }

private:

    HashFunctors<Item> _hash;
    size_t             _nbHash;
    u_int64_t          _size;
    Result*            _blooma;
};

/********************************************************************************/

/* EXPERIMENTAL (not documented). */
template <typename Item, size_t prec=1> class BloomGroup : public system::SmartPointer
{
public:

    class Result
    {
    public:
        Result (u_int64_t v=0)  { memset (v); }

              u_int64_t& operator[] (size_t idx)       { return value[idx]; }
        const u_int64_t& operator[] (size_t idx) const { return value[idx]; }

        const u_int64_t* array () const { return value; }

        Result& operator&= (const Result& r)
        {
            for (size_t j=0; j<prec; j++)  { (*this)[j] &=  r[j]; }
            return *this;
        }

    private:
        u_int64_t value[prec];
        void memset (u_int64_t v)  {  system::impl::System::memory().memset (value, v, prec*sizeof(u_int64_t));  }

        friend class BloomGroup<Item,prec>;
    };

    /** */
    BloomGroup (u_int64_t size, u_int64_t maxMemory, size_t nbHash=4)
        : _hash(nbHash), _nbHash(nbHash), _size(size), _blooma(0)
    {
        printf ("BloomGroup:  size=%ld   sizeof(Result)=%d  maxMemory=%ld\n", size, sizeof(Result), maxMemory);
        if (_size*sizeof(Result) > maxMemory)
        {
            _size = maxMemory /sizeof (Result);
        }
        else
        {
            maxMemory = _size*sizeof(Result);
        }
        printf ("===> size=%ld   allocMemory=%ld\n", _size, sizeof(Result)*_size);

        _blooma = (Result*) MALLOC (_size*sizeof(Result));
        system::impl::System::memory().memset (_blooma, 0, _size*sizeof(Result));
    }

    /** */
    BloomGroup (const std::string& uri)
        : _hash(0), _nbHash(0), _size(0), _blooma(0)
    {
        load (uri);
    }

    /** */
    ~BloomGroup ()  {  FREE (_blooma); }

    /** */
    std::string getName () const { return "BloomGroup"; }

    /** Return the size (in bytes). */
    u_int64_t getMemSize () const { return _size*sizeof(Result); }

    /** */
    void save (const std::string& uri)
    {
        system::IFile* file = system::impl::System::file().newFile (uri, "wb+");
        if (file != 0)
        {
            /** We write the nb of hash functions. */
            file->fwrite (&_nbHash, sizeof(_nbHash), 1);

            /** We write the size of the blooms. */
            file->fwrite (&_size, sizeof(_size), 1);

            /** We write the blooms info. */
            file->fwrite (_blooma, _size*sizeof(Result), 1);

            delete file;
        }
    }

    /** */
    void load (const std::string& uri)
    {
        system::IFile* file = system::impl::System::file().newFile (uri, "rb+");
        if (file != 0)
        {
            /** We read the nb of hash functions. */
            file->fread (&_nbHash, sizeof(_nbHash), 1);

            /** We read the size of the blooms. */
            file->fread (&_size, sizeof(_size), 1);

            /** We allocate the array. */
            _blooma = (Result*) MALLOC (_size*sizeof(Result));
            system::impl::System::memory().memset (_blooma, 0, _size*sizeof(Result));

            /** We read the blooms info. */
            file->fread (_blooma, _size*sizeof(Result), 1);

            delete file;
        }
    }

    /** Insert an item in the 'idx'th Bloom filter
     * \param[in] item : item to be inserted into the filter
     * \param[in] idx : index of the filter */
    void insert (const Item& item, size_t idx)
    {
        u_int64_t q,mask;  euclidian(idx,q,mask);

        for (size_t i=0; i<this->_nbHash; i++)
        {
            u_int64_t h1 = this->_hash (item, i) % this->_size;

#if 1
            this->_blooma[h1][q] |= mask;
#else
            __sync_fetch_and_or (this->_blooma[h1].value + q, mask);
#endif
        }
    }

    /** */
    bool contains (const Item& item, size_t idx)
    {
        u_int64_t q,mask;  euclidian(idx,q,mask);

        for (size_t i=0; i<this->_nbHash; i++)
        {
            u_int64_t h1 = this->_hash (item, i) % this->_size;
            if ( (_blooma[h1][q] & mask) != mask )  {  return false;  }
        }
        return true;
    }

    /** */
    Result contains (const Item& item)
    {
        Result res (~0);

        for (size_t i=0; i<this->_nbHash; i++)
        {
            u_int64_t h1 = this->_hash (item, i) % this->_size;
            res &=  _blooma [h1];
        }
        return res;
    }

private:

    HashFunctors<Item> _hash;
    size_t             _nbHash;
    u_int64_t          _size;
    Result*            _blooma;

    void euclidian (size_t idx, u_int64_t& dividend, u_int64_t& mask) const
    {
        dividend = idx / (8*sizeof(u_int64_t));
        mask     = ((u_int64_t) 1) << (idx % (8*sizeof(u_int64_t)));
    }
};


/********************************************************************************/

/* EXPERIMENTAL (not documented). */
template <typename Item, size_t prec=1> class BloomGroupCacheCoherent : public system::SmartPointer
{
public:

    typedef tools::math::LargeInt<prec> Result;

    /** */
    BloomGroupCacheCoherent (u_int64_t size, size_t nbHash=4, size_t block_nbits=12)
        : _hash(nbHash), _nbHash(nbHash), _size(size), _blooma(0), _nbits_BlockSize(block_nbits)
    {
        _size += (1<<_nbits_BlockSize);

        _blooma = (Result*) MALLOC (_size*sizeof(Result));
        system::impl::System::memory().memset (_blooma, 0, _size*sizeof(Result));

        _mask_block   = (1<<_nbits_BlockSize) - 1;
        _reduced_size = this->_size -  (1<<_nbits_BlockSize) ;
    }

    /** */
    BloomGroupCacheCoherent (const std::string& uri)
        : _hash(0), _nbHash(0), _size(0), _blooma(0)
    {
        load (uri);

        _mask_block   = (1<<_nbits_BlockSize) - 1;
        _reduced_size = this->_size -  (1<<_nbits_BlockSize) ;
    }

    /** */
    ~BloomGroupCacheCoherent ()  {  system::impl::System::memory().free (_blooma); }

    /** */
    std::string getName () const { return "BloomGroupCacheCoherent"; }

    /** Return the size (in bytes). */
    u_int64_t getMemSize () const { return _size*sizeof(Result); }

    /** */
    void save (const std::string& uri)
    {
        system::IFile* file = system::impl::System::file().newFile (uri, "wb+");
        if (file != 0)
        {
            /** We write the nb of hash functions. */
            file->fwrite (&_nbHash, sizeof(_nbHash), 1);

            /** We write the size of the blooms. */
            file->fwrite (&_size, sizeof(_size), 1);

            /** We write block size information. */
            file->fwrite (&_nbits_BlockSize, sizeof(_nbits_BlockSize), 1);

            /** We write the blooms info. */
            file->fwrite (_blooma, _size*sizeof(Result), 1);

            delete file;
        }
    }

    /** */
    void load (const std::string& uri)
    {
        system::IFile* file = system::impl::System::file().newFile (uri, "rb+");
        if (file != 0)
        {
            /** We read the nb of hash functions. */
            file->fread (&_nbHash, sizeof(_nbHash), 1);

            /** We read the size of the blooms. */
            file->fread (&_size, sizeof(_size), 1);

            /** We read block size information. */
            file->fread (&_nbits_BlockSize, sizeof(_nbits_BlockSize), 1);

            /** We allocate the array. */
            _blooma = (Result*) MALLOC (_size*sizeof(Result));
            system::impl::System::memory().memset (_blooma, 0, _size*sizeof(Result));

            /** We read the blooms info. */
            file->fread (_blooma, _size*sizeof(Result), 1);

            delete file;
        }
    }

    /** */
    void insert (const Item& item, size_t idx)
    {
        static const Result ONE (1);

        /** First hash. */
        u_int64_t h0 = this->_hash (item,0) % _reduced_size;
        this->_blooma[h0] |= (ONE << idx);

        for (size_t i=1; i<this->_nbHash; i++)
        {
            /** Other hash. */
            u_int64_t h1 = h0  + (simplehash16 (item,i) & _mask_block);
            this->_blooma[h1] |= (ONE << idx);
        }
    }

    /** */
    bool contains (const Item& item, size_t idx)
    {
        static const Result ONE  (1);

        /** First hash. */
        u_int64_t h0 = this->_hash (item,0) % _reduced_size;
        if ((this->_blooma[h0] & (ONE << idx)) != (ONE << idx))  {  return false;  }

        for (size_t i=1; i<this->_nbHash; i++)
        {
            /** Other hash. */
            u_int64_t h1 = h0  + (simplehash16 (item,i) & _mask_block);
            if ((this->_blooma[h1] & (ONE << idx)) != (ONE << idx))  {  return false;  }
        }
        return true;
    }

    /** */
    Result contains (const Item& item)
    {
        static const Result ZERO (0);
        Result res = ~ZERO;

        u_int64_t h0 = this->_hash (item,0) % _reduced_size;
        res &= _blooma [h0];

        for (size_t i=1; i<this->_nbHash; i++)
        {
            u_int64_t h1 = h0 + (simplehash16(item,i) & _mask_block);
            res &= _blooma [h1];
        }
        return res;
    }

private:

    HashFunctors<Item> _hash;
    size_t             _nbHash;
    u_int64_t          _size;
    Result*            _blooma;

    u_int64_t  _mask_block;
    size_t     _nbits_BlockSize;
    u_int64_t  _reduced_size;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BLOOM_GROUP_HPP_ */
