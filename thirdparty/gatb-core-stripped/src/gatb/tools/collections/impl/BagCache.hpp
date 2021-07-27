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

/** \file BagCache.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Bag implementation as a cache for a referred Bag
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BAG_CACHE_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BAG_CACHE_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Bag.hpp>
#include <gatb/system/impl/System.hpp>
#include <vector>
#include <algorithm>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
/** Implementation of API from collections */
namespace impl          {
/********************************************************************************/

/** \brief Bag implementation as a cache to a referred Bag instance
 *
 * This implementation is a proxy to a Bag instance. One can insert items into such
 * a bag, it will just put them into a cache. When the cache is full, all the cached
 * items are inserted into the delegate; this operation may be protected by a
 * synchronizer in case several threads use different BagCache on the same delegate.
 */
template <typename Item> class BagCache : public Bag<Item>, public system::SmartPointer
{
public:

    /** Constructor */
    BagCache () : _ref(0), _synchro(0), _items(0), _nbMax(0), _idx(0)  {}

    /** Constructor. Cache size is in terms of number of items (I think) */
    BagCache (Bag<Item>* ref, size_t cacheSize, system::ISynchronizer* synchro=0)
        : _ref(0), _synchro(0), _items(0), _nbMax(cacheSize), _idx(0)
    {
        setRef     (ref);
        setSynchro (synchro);
        _items = (Item*) CALLOC (_nbMax, sizeof(Item));
        system::impl::System::memory().memset (_items, 0, _nbMax*sizeof(Item));
    }

    BagCache (const BagCache<Item>& b)
        : _ref(0), _synchro(0), _items (0), _nbMax(b._nbMax), _idx(0)
    {
        setRef     (b._ref);
        setSynchro (b._synchro);

        _items = (Item*) CALLOC (_nbMax, sizeof(Item));
        system::impl::System::memory().memset (_items, 0, _nbMax*sizeof(Item));
    }

    /** Destructor. */
    virtual ~BagCache ()
    {
        /** We flush the potential remaining stuf. */
        flush ();

        /** We clean resources. */
        setRef     (0);
        setSynchro (0);

        FREE (_items);
    }

    /**  \copydoc Bag::insert */
    void insert (const Item& item)
    {
        if (_idx+1 > _nbMax)
        {
            if (_synchro)  {  _synchro->lock();    }
            flushCache ();
            if (_synchro)  {  _synchro->unlock();  }
        }

        _items[_idx++] = item;
    }

    /**  \copydoc Bag::flush */
    void flush ()
    {
        if (_synchro)  {  _synchro->lock();    }
        flushCache ();
        _ref->flush();
        if (_synchro)  {  _synchro->unlock();  }
    }

protected:

    Bag<Item>*             _ref;
    void setRef (Bag<Item>* ref)  { SP_SETATTR(ref); }

    system::ISynchronizer* _synchro;
    void setSynchro (system::ISynchronizer* synchro)  { SP_SETATTR(synchro); }

    Item*                  _items;
    size_t                 _nbMax;
    size_t                 _idx;

    void flushCache ()
    {
        if (_idx > 0)  { _ref->insert (_items, _idx);  }
        _idx = 0;
    }
};
    
/********************************************************************************/

/** \brief Bag implementation as a cache to a referred Bag instance.
 *
 * The cache is sorted before sent to the reference
 */
template <typename Item> class BagCacheSortedBuffered : public BagCache<Item>
{
public:
    
    /** Constructor. */
    BagCacheSortedBuffered (Bag<Item>* ref, size_t cacheSize, Item* sharedBuffer, size_t sharedCacheSize,  size_t * idxShared,  system::ISynchronizer* outsynchro=0, system::ISynchronizer* synchro=0
)
    : BagCache<Item>(ref, cacheSize, synchro), _outsynchro(0),_idxShared(idxShared)  {
        _sharedBuffer = sharedBuffer;
        setOutSynchro (outsynchro);
        _sharedCacheSize = sharedCacheSize;
      //  printf("BagCacheSortedBuffered created with _idxShared %p : %zu \n",_idxShared,*_idxShared);
    }

    
    /** Destructor. */
    virtual ~BagCacheSortedBuffered ()
    {
        /** We flush the potential remaining stuf. */
        //flush ();
        
        /** We clean resources. */
        setOutSynchro (0);
    }
    
    /**  \copydoc Bag::insert */
    void insert (const Item& item)
    {
        if (this->_idx+1 > this->_nbMax)
        {
            if (this->_synchro)  {  this->_synchro->lock();    }
            flushLocalCache ();
            if (this->_synchro)  {  this->_synchro->unlock();  }
        }
        this->_items[this->_idx++] = item;
    }
    
    /**  \copydoc Bag::flush */
    void flush ()
    {
        //printf("main flush \n");
        if (this->_synchro)  {  this->_synchro->lock();    }
        
        flushLocalCache ();

        flushCache ();
        if (this->_synchro)  {  this->_synchro->unlock();  }

        if (this->_outsynchro)  {  this->_outsynchro->lock();    }
        this->_ref->flush();
        if (this->_outsynchro)  {  this->_outsynchro->unlock();  }
    }
    
private:
    Item * _sharedBuffer; // shared buffer for sorting, allocated from outside
    system::ISynchronizer* _outsynchro;
    void setOutSynchro (system::ISynchronizer* outsynchro)  { SP_SETATTR(outsynchro); }
    size_t _sharedCacheSize;
    size_t * _idxShared;

    void flushLocalCache () //flush local cache to shared buffer
    {
      //  printf("flush local cache to %p  from %p  %zu \n",_sharedBuffer + *_idxShared,this->_items,this->_idx);
        if (this->_idx > 0)  {
            
            while(*_idxShared + this->_idx > _sharedCacheSize)
            // must be a while, not an 'if' : flushCache may release temporarily the lock, so must re-check when exiting flushcache
            {
                //sort and flush shared buffer to ref bag
                flushCache ();
            }

            system::impl::System::memory().memcpy (_sharedBuffer + *_idxShared, this->_items, this->_idx*sizeof(Item));
            *_idxShared += this->_idx;
            
            //this->_ref->insert (this->_items, this->_idx);

        }
        
        this->_idx = 0;
    }
    
    void flushCache () //flush shared buffer to ref bag
    {
     //   printf("flush  shared cache %zu  from  %p this %p\n",*_idxShared,_sharedBuffer,this);

        if (*_idxShared > 0)
        {
            //first duplicate into temp vector
            std::vector<Item> tempArray(_sharedBuffer, _sharedBuffer + *_idxShared);
            *_idxShared = 0;
            //here shared buffer is effectively flushed to temp array, safe to unlock

            if (this->_synchro)  {  this->_synchro->unlock();    }
            std::sort (tempArray.begin(), tempArray.end());

            //flush to file using outsynchro
            if (this->_outsynchro)  {  this->_outsynchro->lock();    }
            this->_ref->insert (tempArray, tempArray.size());
            if (this->_outsynchro)  {  this->_outsynchro->unlock();  }

            if (this->_synchro)  {  this->_synchro->lock();    }
        }
        //  if (*_idxShared > 0)  { this->_ref->insert (_sharedBuffer, *_idxShared);  }
        //  *_idxShared = 0;
    }
};
    
/********************************************************************************/

/** \brief Bag implementation as a cache to a referred Bag instance.
 *
 * The cache is sorted before sent to the reference
 */
template <typename Item> class BagCacheSorted : public BagCache<Item>
{
public:

    /** Constructor. */
    BagCacheSorted (Bag<Item>* ref, size_t cacheSize, system::ISynchronizer* synchro=0)
        : BagCache<Item>(ref, cacheSize, synchro)  { }

    /**  \copydoc Bag::insert */
    void insert (const Item& item)
    {
        if (this->_idx+1 > this->_nbMax)
        {
            std::sort (this->_items.begin (), this->_items.end ());

            if (this->_synchro)  {  this->_synchro->lock();      }
            this->flushCache ();
            if (this->_synchro)  {  this->_synchro->unlock();    }
        }

        this->_items[this->_idx++] = item;
    }

    /**  \copydoc Bag::flush */
    void flush ()
    {
        this->_items.resize (this->_idx);
        std::sort (this->_items.begin (), this->_items.end ());

        if (this->_synchro)  {  this->_synchro->lock();    }
        this->flushCache ();
        this->_ref->flush();
        if (this->_synchro)  {  this->_synchro->unlock();  }
    }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_BAG_CACHE_HPP_ */
