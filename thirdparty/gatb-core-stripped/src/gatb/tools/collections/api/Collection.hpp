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

/** \file Collection.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Collection interface
 *
 *  This file holds interfaces related to the Collection interface
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_COLLECTION_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_COLLECTION_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/collections/api/Bag.hpp>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
/********************************************************************************/

/** \brief Collection interface
 *
 * The Collection interface is the union of a Bag and an Iterable interfaces
 *
 * It is also to get/set properties (as [key,value]) to collections.
 */
template <class Item> class Collection : public Bag<Item>, public Iterable<Item>
{
public:

    /** Destructor. */
    virtual ~Collection () {}

    /** \return the bag instance. */
    virtual Bag<Item>* bag() = 0;

    /** \return the iterable instance. */
    virtual Iterable<Item>* iterable() = 0;

    /** Remove physically the collection. */
    virtual void remove () = 0;

    /** Add a property to the collection.
     * \param[in] key : key of the property
     * \param[in] value  : value of the property. */
    virtual void addProperty (const std::string& key, const std::string value) = 0;

    /** Add a property to the collection.
     * \param[in] key : key of the property
     * \param[in] format  : printf like arguments
     * \param[in] ... : variable arguments. */
    virtual void addProperty (const std::string& key, const char* format, ...) = 0;

    /** Retrieve a property for a given key
     * \param[in] key : key of the property to be retrieved
     * \return the value of the property. */
    virtual std::string getProperty (const std::string& key) = 0;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_COLLECTION_HPP_ */
