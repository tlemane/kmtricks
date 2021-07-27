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

/** \file Bag.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Bag interface
 *
 *  This file holds interfaces related to the Bag interface, ie something we can put items into.
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_BAG_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_BAG_HPP_

/********************************************************************************/

#include <gatb/system/api/ISmartPointer.hpp>
#include <gatb/tools/designpattern/api/Iterator.hpp>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
/********************************************************************************/

/** \brief Bag interface
 *
 * A bag is something we can put items into. Such items are typed and the type is
 * provided as a template.
 *
 * Several ways for inserting items exist. There is also a flush method that makes
 * sure that all the inserted items are actually in the bag.
 */
template <class Item> class Bag : public virtual  system::ISmartPointer
{
public:

    /** Insert an item into the bag.
     * \param[in] item : the item to be inserted. */
    virtual void insert (const Item& item) = 0;

    /** Insert items into the bag.
     * \param[in] items : items to be inserted.
     * \param[in] length : the number of items to be inserted. If 0 (default value),
     *                     all the items of the vector are inserted. */
    virtual void insert (const std::vector<Item>& items, size_t length=0)
    {
        size_t n = length==0 ? items.size() : length;
        for (size_t i=0; i<n; i++)  {  insert (items[i]); }
    }

    /** Insert items into the bag.
     * \param[in] items : items to be inserted.
     * \param[in] length : number of items to be inserted. */
    virtual void insert (const Item* items, size_t length)
    {
        for (size_t i=0; i<length; i++)  {  insert (items[i]); }
    }

    /** Flush the current content. May be useful for implementation that uses a cache. */
    virtual void flush () = 0;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_BAG_HPP_ */
