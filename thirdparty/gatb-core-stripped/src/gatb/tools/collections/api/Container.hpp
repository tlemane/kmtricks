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

/** \file Container.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Container interface
 *
 *  This file holds interfaces related to the Container interface, ie something we can ask for items.
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_CONTAINER_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_CONTAINER_HPP_

/********************************************************************************/

#include <gatb/system/api/ISmartPointer.hpp>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
/********************************************************************************/

/** \brief Container interface
 *
 * The Container interface provides an operation that ask for a given item
 */
template <class Item> class Container : public virtual system::ISmartPointer
{
public:

    /** Destructor. */
    virtual ~Container() {}

    /** Tells whether an item exists or not
     * \return true if the item exists, false otherwise */
    virtual bool contains (const Item& item) = 0;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_CONTAINER_HPP_ */
