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

/** \file Iterable.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Iterable interface
 *
 *  This file holds interfaces related to the Iterable interface
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_ITERABLE_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_ITERABLE_HPP_

/********************************************************************************/

#include <gatb/system/api/ISmartPointer.hpp>
#include <gatb/tools/designpattern/api/Iterator.hpp>

/********************************************************************************/
namespace gatb          {
namespace core          {
/** \brief Tools package */
namespace tools         {
/** \brief Collections interfaces */
namespace collections   {
/********************************************************************************/

/** \brief Iterable interface
 *
 * The Iterable interface provides an operation that creates an iterator. It also
 * provides other methods that give the exact number (or an estimation) of items.
 *
 * Note that one Iterable instance can create several iterators.
 */
template <class Item> class Iterable : public virtual system::ISmartPointer
{
public:

    virtual ~Iterable() {}

    /** Create an iterator for the given Iterable instance.
     * \return the new iterator. */
    virtual dp::Iterator<Item>* iterator () = 0;

    /** Direct iteration through a functor.
     * \param[in] f : the functor to be applied on each sequence of the bank. */
    template<typename Functor> void iterate (Functor f)
    {
        tools::dp::Iterator<Item>* it = this->iterator();  LOCAL (it);
        for (it->first(); !it->isDone(); it->next())  { f((*it).item()); }
    }

    /** Return the number of items. If a specific implementation doesn't know the value,
     * it should return -1 by convention.
     * \return the number of items if known, -1 otherwise. */
    virtual int64_t getNbItems () = 0;

    /** Return the (approximate) number of items. If a specific implementation doesn't know the value,
     * it should return -1 by convention.
     * \return the number of items if known, -1 otherwise. */
    virtual int64_t estimateNbItems () = 0;

    /** Return a buffer of items.
     * \param[out] buffer : the buffer
     * \return the buffer */
    virtual Item* getItems (Item*& buffer)
    {
        std::cout << "ERROR" << std::endl; exit(1);
        //throw "Iterable::getItems... SHOULD NOT BE HERE..."; // this throw wasn't caught by the unit test. so i'm doing an exit now.
        return buffer;
    }

    /** Return a buffer of items.
     * \param[out] buffer : the buffer
     * \param[in] start : index where to start in the buffer
     * \param[in] nb : number of items to be retrieved
     * \return the number of items retrieved */
    virtual size_t getItems (Item*& buffer, size_t start, size_t nb)
    {
        std::cout << "ERROR" << std::endl; exit(1);
        //throw "Iterable::getItems... SHOULD NOT BE HERE...";
        return 0;
    }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_ITERABLE_HPP_ */
