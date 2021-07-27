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

/** \file CollectionAbstract.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Collection interface
 *
 *  This file holds interfaces related to the Collection interface
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_COLLECTION_ABSTRACT_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_COLLECTION_ABSTRACT_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Collection.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>
#include <gatb/system/impl/System.hpp>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

/** \brief Abstract implementation of the Collection interface.
 *
 * This class implements the Collection interface by delegating the job to an instance
 * of Bag and an instance of Iterable.
 *
 * All the methods are delegated to one of these two instances.
 */
template <class Item>
class CollectionAbstract : public Collection<Item>
{
public:

    /** Constructor.
     * \param bag      : reference on the bag delegate.
     * \param iterable : reference on the iterable delegate
     */
    CollectionAbstract (Bag<Item>* bag, Iterable<Item>* iterable)
        : _bag(0), _iterable(0)
    {
        setBag      (bag);
        setIterable (iterable);
    }

    /** Destructor. */
    virtual ~CollectionAbstract()
    {
        setBag      (0);
        setIterable (0);
    }

    /** \copydoc Collection::bag */
    Bag<Item>* bag() { return _bag; }

    /** \copydoc Collection::iterable */
    Iterable<Item>* iterable()  { return _iterable; }

    /** \copydoc Iterable::iterator */
    dp::Iterator<Item>* iterator ()  { return _iterable->iterator(); }

    /** \copydoc Iterable::getNbItems */
    int64_t getNbItems ()  { return _iterable->getNbItems(); }

    /** \copydoc Iterable::estimateNbItems */
    int64_t estimateNbItems () { return _iterable->estimateNbItems(); }

    /** \copydoc Iterable::getItems */
    Item* getItems (Item*& buffer)  { 
        //std::cout << "CollectionAsbtract getItems called" << std::endl;
        return _iterable->getItems(buffer); }

    /** \copydoc Iterable::getItems(Item*& buffer, size_t start, size_t nb) */
    size_t getItems (Item*& buffer, size_t start, size_t nb)  { 
        //std::cout << "CollectionAsbtract getItems called" << std::endl;
        return _iterable->getItems (buffer, start, nb); }

    /** \copydoc Bag::insert */
    void insert (const Item& item)  { _bag->insert (item); }

    /** \copydoc Bag::insert(const std::vector<Item>& items, size_t length) */
    void insert (const std::vector<Item>& items, size_t length)  { _bag->insert (items, length); }

    /** \copydoc Bag::insert(const Item* items, size_t length) */
    void insert (const Item* items, size_t length)  { _bag->insert (items, length); }

    /** \copydoc Bag::flush */
    void flush ()  { _bag->flush(); }

    /** \copydoc Collection::addProperty */
    void addProperty (const std::string& key, const std::string value)  {
        std::cout << "warning: collectionAbstract.addProperty() called without an implementation (notify a developer)" << std::endl;
    }

    /** \copydoc Collection::addProperty */
    void addProperty (const std::string& key, const char* format ...)
    {
        va_list args;
        va_start (args, format);
        this->addProperty (key, misc::impl::Stringify::format(format, args));
        va_end (args);
    }

    /** \copydoc Collection::getProperty */
    std::string getProperty (const std::string& key)  {  return std::string("");  }

protected:

    Bag<Item>* _bag;
    void setBag (Bag<Item>* bag)  { SP_SETATTR(bag); }

    Iterable<Item>* _iterable;
    void setIterable (Iterable<Item>* iterable)  { SP_SETATTR(iterable); }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_COLLECTION_ABSTRACT_HPP_ */
