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

/** \file CollectionCache.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Collection interface
 *
 *  This file holds interfaces related to the Collection interface
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_COLLECTION_CACHE_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_COLLECTION_CACHE_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Collection.hpp>
#include <gatb/tools/collections/impl/BagCache.hpp>
#include <gatb/system/impl/System.hpp>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

/** \brief Cache implementation of the Collection interface
 *
 * This implementations allows to add items into the collection through a cache. It
 * means that a memory block for inserting the items and when the memory is full,
 * all its items are flushed into the delegate Collection (provided at construction).
 *
 * A synchronizer can be provided to be sure that the delegate can't suffer concurrent
 * accesses by several threads. This allows to encapsulate a single Collection into several
 * CollectionCache instances, each one used in a thread.
 */
template <class Item> class CollectionCache : public CollectionAbstract<Item>, public system::SmartPointer
{
public:

    /** Constructor. */
    CollectionCache (Collection<Item>& ref,  size_t cacheSize, system::ISynchronizer* synchro)
        : CollectionAbstract<Item> (
            new BagCache<Item> (ref.bag(), cacheSize, synchro),
            ref.iterable()
        ), _ref(ref)  {}

    /** Destructor. */
    virtual ~CollectionCache() {}

    /** \copydoc Collection::remove */
    void remove ()  { _ref.remove(); }

    /** Accessor to the delegate Collection.
     * \return the delegate Collection instance. */
    Collection<Item>& getRef ()  { return _ref; }

private:
    Collection<Item>& _ref;
};

/********************************************************************************/
    
/** \brief Cache implementation of the Collection interface
 *
 * This implementations is like CollectionCache. The difference is that
 * the items are first sorted before inserting them into the delegate
 * Collection instance.
 *
 * Remark : maybe could be inherited from CollectionCache ?
 */
template <class Item> class CollectionCacheSorted : public CollectionAbstract<Item>, public system::SmartPointer
{
public:
    
    /** Constructor. */
    CollectionCacheSorted (Collection<Item>& ref,  size_t cacheSize, size_t sharedCacheSize,  system::ISynchronizer* synchro, system::ISynchronizer* outsynchro, Item* sharedBuffer, size_t * idxShared) //
    : CollectionAbstract<Item> (
            new BagCacheSortedBuffered<Item> (
                ref.bag(),
                cacheSize,
                sharedBuffer,
                sharedCacheSize,
                idxShared,
                outsynchro,synchro
            ),
            ref.iterable()
        ),
        _ref(ref)  {}
    
    /** Destructor. */
    virtual ~CollectionCacheSorted() {}
    
    /** \copydoc Collection::remove */
    void remove ()  { _ref.remove(); }
    
    /** Accessor to the delegate Collection.
     * \return the delegate Collection instance. */
    Collection<Item>& getRef ()  { return _ref; }
    
private:
    Collection<Item>& _ref;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_COLLECTION_CACHE_HPP_ */
