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

/** \file Observer.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Implementation of INode interface
 */

#ifndef _GATB_CORE_DP_ITERATOR_IMPL_OBSERVER_HPP_
#define _GATB_CORE_DP_ITERATOR_IMPL_OBSERVER_HPP_

/********************************************************************************/

#include <gatb/tools/designpattern/api/IObserver.hpp>
#include <list>

/********************************************************************************/
namespace gatb  {
namespace core  {
namespace tools {
namespace dp    {
namespace impl  {
/********************************************************************************/


/** \brief Class that notifies potential observers.
 *
 * The main purpose of this class is to manage the set of IObservers attached to the subject.
 *
 * Then, classes that want subject-like behavior can inherit from Subject or have a Subject
 * attribute.
 *
 * \see ISubject
 */
class Subject : public ISubject
{
public:

    /** Constructor. */
    Subject ();

    /** Constructor.
     * \param[in] interface : the identifier of the subject. */
    Subject (const InterfaceId& interface);

    /** Destructor. */
    virtual ~Subject();

    /** Returns the identifier of the subject. */
    InterfaceId getInterface ()  { return _interface; }

    /** Attach an observer to this subject.
     * \param[in] observer : the observer to be attached.
     */
    void addObserver    (IObserver* observer);

    /** Detach an observer from this subject.
     * \param[in] observer : the observer to be detached.
     */
    void removeObserver (IObserver* observer);

    /** Notify observers by sending a EventInfo instance.
     * \param[in]  event : the information to be sent to the observers.
     */
    void notify (EventInfo* event);

private:

    /** Identifier of the subject. */
    InterfaceId             _interface;

    /** Optimization attribute. */
    IObserver*              _singleObserver;

    /** List of observers attached to the subject. */
    std::list<IObserver*>*  _observers;
};


/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DP_ITERATOR_IMPL_CELL_HPP_ */
