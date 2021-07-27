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

/** \file IterableHelpers.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Helpers class for iterable concept
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_ITERABLE_HELPERS_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_ITERABLE_HELPERS_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

/** \brief Iterator with progress information.
 *
 * This Iterator implementation encapsulates a given iterator and shows progress
 * information to end users.
 *
 * The way the progress information is displayed is set through the second template
 * argument of the class which is some IteratorListener implementation class. There
 * is a default argument for this template, so one has often have only to give the
 * first template type (ie the kind of items to be iterated).
 *
 * This is an example of the Decorator design pattern.
 */
template<class Type, class Listener=tools::misc::impl::ProgressDefault>
class ProgressIterator : public tools::dp::impl::SubjectIterator<Type>
{
public:

    /** Constructor. It uses the provided Iterable object to build the delegate Iterator
     * instance.
     * \param[in] iterable : used to create the delegate iterator
     * \param[in] msg : message displayed at each progression notification
     * \param[in] divide : number of notifications to be send
     */
    ProgressIterator (Iterable<Type>& iterable, const char* msg = "progress", size_t divide=100)
        : tools::dp::impl::SubjectIterator<Type> (
            iterable.iterator(),
            (iterable.getNbItems() >= 0 ? iterable.getNbItems() : iterable.estimateNbItems()) / divide,
            new Listener ((iterable.getNbItems() >= 0 ? iterable.getNbItems() : iterable.estimateNbItems()), msg)
    ) {}

    /** Constructor. It uses the provided Iterator object as the delegate iterator
     * \param[in] iterator : delegate iterator
     * \param[in] msg : message displayed at each progression notification
     * \param[in] nbItems : number of items to be iterated
     */
    ProgressIterator (tools::dp::Iterator<Type>* iterator, const char* msg, size_t nbItems)
        : tools::dp::impl::SubjectIterator<Type> (
            iterator,
            nbItems / 100,
            new Listener (nbItems, msg)
    ) {}
};

/********************************************************************************/

/** \brief Adaptor of an Iterable<T1> into an Iterable<T2>
 *
 * This class converts the type T1 of an iterable into another iterable of type T2.
 *
 * A functor must be provided in order to convert one item of T1 into T2.
 */
template <class T1, class T2, class Adaptor>
class IterableAdaptor : public Iterable<T2>, public system::SmartPointer
{
public:
    /** */
    IterableAdaptor (Iterable<T1>& ref)  : _ref(ref) {}

    /** Create an iterator for the given Iterable instance.
     * \return the new iterator. */
    dp::Iterator<T2>* iterator ()  { return new tools::dp::impl::IteratorAdaptor<T1,T2,Adaptor> (_ref.iterator()); }

    /** Return the number of items. If a specific implementation doesn't know the value,
     * it should return -1 by convention.
     * \return the number of items if known, -1 otherwise. */
    int64_t getNbItems () { return _ref.getNbItems(); }

    /** Return the (approximate) number of items. If a specific implementation doesn't know the value,
     * it should return -1 by convention.
     * \return the number of items if known, -1 otherwise. */
    int64_t estimateNbItems ()  { return _ref.estimateNbItems(); }

    /** Return a buffer of items.
     * \param[out] buffer : the buffer
     * \return the buffer */
    T2* getItems (T2*& buffer) { std::cout << "IterableAdaptor::getItems() called (but it's not implemented)" << std::endl; exit(1); return 0; }

    /** */
    size_t getItems (T2*& buffer, size_t start, size_t nb) { std::cout << "IterableAdaptor::getItems() called (but that function isn't implemented)" << std::endl; exit(1); return 0; }

private:
    Iterable<T1>& _ref;
};

/********************************************************************************/

/** \brief some helper methods for Iterable objects. */
class IterableHelpers
{
public:

    /** Retrieve in a vector items from an Iterable instance.
     * \param[in] iterable : the instance we want to get items from
     * \param[out] items : vector of items to be retrieved. This vector must be sized to
     * the number of items that one wants to retrieve. */
    template<typename T>
    static bool getItems (Iterable<T>& iterable, std::vector<T>& items)
    {
        std::cout << "IterableHelper::getItems() called" << std::endl; 
        size_t nbItems = items.size();

        dp::Iterator<T>* it = iterable.iterator();  LOCAL (it);
        it->get (items);

        return items.size() == nbItems;
    }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_ITERABLE_HELPERS_HPP_ */
