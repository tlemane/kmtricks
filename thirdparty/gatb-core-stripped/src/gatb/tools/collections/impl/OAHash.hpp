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

/** \file OAHash.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Container implementation
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_OAHASH_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_OAHASH_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Container.hpp>
#include <gatb/tools/collections/api/Bag.hpp>
#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/api/Abundance.hpp>
#include <algorithm>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

/** \brief Hash table implementation
 */
template <typename Item> class OAHash
{
    typedef misc::Abundance<Item> element_pair;

public:

    /** Get the size (in byte) of an item.
     * \return the item size.
     */
    static int size_entry ()  {  return sizeof(element_pair); }

    /** Get the max number of items for the hash table
     * \return the max items number.
     */
    int getMaxNbItems ()  { return hash_size; }

    /** Constructor.
     * \param[in] max_memory : max memory for the hash table.
     */
    OAHash (u_int64_t max_memory)
    {
        hash_size = max_memory / sizeof(element_pair);
        if (hash_size == 0)  {  throw system::Exception ("empty OAHash allocated");  }
        data = (element_pair *) CALLOC ( hash_size, sizeof(element_pair));  //create hashtable
    }

    /** Destructor. */
    ~OAHash()  {  FREE (data);  }

    /** Insert an item with its value into the hash table.
     * \param[in] graine : key
     * \param[in] value : value
     */
    void insert (const Item& graine, int value)
    {
        element_pair *element = find_slot(graine);
        if (!is_occupied(element))
            element->value = graine;
        element->abundance = value;
    }

    /** Increment the value for a given key.
     * \param[in] graine : key
     */
    void increment (const Item& graine)
    {
        element_pair *element = find_slot(graine);

        if (!is_occupied(element))
            element->value = graine;
        element->abundance = element->abundance + 1;
    }

    /** Get the value for a given key
     * \param[in] graine : key
     * \param[out] val : value to be retrieved
     * \return true if the key is found, false otherwise.
     */
    bool get (const Item& graine, int * val=0)
    {
        element_pair *element = find_slot(graine, false);
        if (element == 0)  { return 0; }

        if (!is_occupied(element))
            return false;
        if (element->value == graine)
            if (val != NULL)
                *val = element->abundance;
        return true;
    }

    /** Tells whether or not the given key is in the table
     * \param[in] graine : key to be checked
     * \return true if present, false otherwise.
     */
    bool has_key (const Item& graine)  {     return get(graine,NULL) == 1;  }

    /** Get the memory usage of the hash table
     * \return the memory usage (in bits).
     */
    u_int64_t memory_usage() {     return hash_size* sizeof(element_pair); /* in bits */ }

    /** Get the load factor of the hash table
     * \return the load factor.
     */
    float load_factor()
    {
        u_int64_t ptr = 0;
        u_int64_t nbKeys = 0;
        while (ptr < hash_size)
        {
            while ((ptr < hash_size) &&  ((data+ptr)->abundance == 0)  )
                ptr++;

            if (ptr == hash_size)
                break;

            nbKeys++;
            ptr++;
        }
        return (float)nbKeys/(float)hash_size;
    }

    /** Get an iterator for the hash table.
     * \param[in] sorted : if true, items are iterated in a sorted way
     * \return an iterator over the items of the hash table.
     */
    dp::Iterator <misc::Abundance<Item> >* iterator (bool sorted=false)
    {  if (sorted==false) { return new Iterator(*this); } else { return new IteratorSorted (*this);  } }


    /************************************************************/
    class Iterator : public tools::dp::Iterator <misc::Abundance<Item> >
    {
    public:

        Iterator (OAHash<Item>& aRef) : ref(aRef), iterator(0), iteratorMax(0), done(true)  {}

        /** \copydoc tools::dp::Iterator::first */
        void first()
        {
            iterator    = ref.data - 1;
            iteratorMax = ref.data + ref.hash_size;
            done        = false;

            next ();
        }

        /** \copydoc tools::dp::Iterator::next */
        void next()
        {
            while (!done)
            {
                ++iterator;
                done = (iterator >= iteratorMax);
                if (!done && iterator->abundance != 0)
                {
                    *this->_item = *iterator;  break;
                }
            }
        }

        /** \copydoc tools::dp::Iterator::isDone */
        bool isDone ()   {  return done; }

        /** \copydoc tools::dp::Iterator::item */
        misc::Abundance<Item>& item ()     { return *this->_item; }

    private:
        OAHash<Item>&  ref;
        element_pair*  iterator;
        element_pair*  iteratorMax;
        bool           done;
    };

    /************************************************************/
    class IteratorSorted : public tools::dp::Iterator <misc::Abundance<Item> >
    {
    public:

        IteratorSorted (OAHash<Item>& aRef) : _ref(aRef), _idx(0), _nb(0)
        {
            if (_ref.hash_size > (1ULL<<32))  { throw system::Exception ("OAHash::sort  too many items..."); }

            for (u_int64_t idx=0; idx<_ref.hash_size; idx++)
            {
                if (_ref.data[idx].abundance != 0)  { _offsets.push_back (idx); }
            }

            _nb = _offsets.size();

            std::sort (_offsets.begin(), _offsets.end(), CmpPair(_ref,_offsets));
        }

        /** \copydoc tools::dp::Iterator::first */
        void first()  {  _idx = -1;  next ();  }

        /** \copydoc tools::dp::Iterator::next */
        void next()  {  ++_idx;  if (_idx < _nb ) { *(this->_item) = (_ref.data[_offsets[_idx]]); }
        }

        /** \copydoc tools::dp::Iterator::isDone */
        bool isDone ()   {  return _idx >= _nb; }

        /** \copydoc tools::dp::Iterator::item */
        misc::Abundance<Item>& item ()     { return *this->_item; }

    private:
        OAHash<Item>&           _ref;
        std::vector<u_int32_t>  _offsets;
        int64_t                 _idx;
        int64_t                 _nb;

        struct CmpPair
        {
            element_pair*           _data;
            std::vector<u_int32_t>& _offsets;
            CmpPair (OAHash<Item>& ref, std::vector<u_int32_t>& offsets) : _data(ref.data), _offsets(offsets) {}
            bool operator() (u_int32_t i1, u_int32_t i2)  {  return _data[i1].value < _data[i2].value;  }
        };
    };

protected:

    u_int64_t      hash_size;
    element_pair*  data;

    /** */
    element_pair * find_slot (const Item& key, bool exceptionOnBadKey = true)
    {
        u_int64_t ptr = oahash (key) % hash_size;
        element_pair* element = data+ptr;
        u_int64_t retries = 0;

        // search until we either find the key, or find an empty slot.
        while ( ( is_occupied(element)) && ( element->value != key ) && (retries < hash_size))
        {
            ptr = (ptr + 1) % hash_size;
            element = data+ptr;
            retries++;
        }

        if (retries == hash_size)
        {
            if (exceptionOnBadKey)  {  throw system::Exception ("OAHash: max rehashes reached: %lld (notify a developer)", hash_size);  }
            return 0;
        }

        return element;
    }


    bool is_occupied (element_pair *element)   {  return (element->abundance != 0); }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_OAHASH_HPP_ */
