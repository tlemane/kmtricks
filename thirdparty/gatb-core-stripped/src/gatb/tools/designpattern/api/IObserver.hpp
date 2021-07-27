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

/** \file IObserver.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Smart Pointer Design Pattern definition
 *
 *  This file holds interfaces related to the Design Pattern Observer.
 */

#ifndef _GATB_CORE_TOOLS_DP_IOBSERVER_HPP_
#define _GATB_CORE_TOOLS_DP_IOBSERVER_HPP_

/********************************************************************************/

#include <gatb/system/api/ISmartPointer.hpp>
#include <gatb/system/api/types.hpp>

/********************************************************************************/
namespace gatb  {
namespace core  {
/** \brief Tools package */
namespace tools {
/** \brief Design Patterns concepts */
namespace dp    {
/********************************************************************************/

/* Forward declarations. */
class ISubject;

/********************************************************************************/

/* We define a type for identifying an interface. */
typedef u_int32_t InterfaceId;

/********************************************************************************/

/** \brief Piece of information linked to a notification.
 *
 *  The EventInfo class is intended to be received by clients during a notification.
 *  Subclasses are meant to provide specific information related to specific notifications.
 *
 *  Defined as a SmartPointer for easing the EventInfo instances life cycle management.
 */
class EventInfo : public system::SmartPointer
{
public:

    /** Returns the identifier of the group the EventInfo belongs to.
     *  \return the identifier
     */
    InterfaceId getInterface ()  const { return _interface; }

protected:

    /** Constructor.
     * \param[in] interface : identifier of the interface this EventInfo instance belongs to.
     */
    EventInfo (const InterfaceId& interface) : _interface(interface) {}

    /** Destructor. */
    virtual ~EventInfo() {}

    /** Id of the group the EventInfo is related. This allows to group several EventInfo
     *  subclasses with a common identifier. */
    InterfaceId _interface;
};

/********************************************************************************/

/** \brief Interface of the Observer Design Pattern.
 *
 *  The Observer Design Pattern is an object oriented way for managing a notification
 *  system in a software. This is an extremely useful pattern because it allows
 *  the communication between software components without having strong dependencies
 *  between them.
 *
 *  For instance, a tool sends notifications about its progression. This is useful for
 *  potential clients that want to know whether they have to wait a long time before
 *  the algorithm finishes. The important point is that a tool is not polluted by clients
 *  concerns, it just notifies some information, without even know who in particular will
 *  receive the information.
 *
 *  This interface merely defines an 'update' method that can be called by some subjects
 *  for telling that some notification is occurring. Some information can be provided to
 *  the IObserver instance through an EventInfo instance.
 *
 *  Note that the subject that initiated the notification is also provided to the observer.
 *  This may be useful for discriminating some unwanted behavior, for instance the IObserver
 *  being itself the subject.
 *
 *  \see ISubject
 *  \see EventInfo
 *  \see IterationStatusEvent
 */
class IObserver
{
public:

    /** Destructor. */
    virtual ~IObserver() {}

    /** Method called by Subject for some notification.
     *  \param[in] evt     : specific information related to the notification.
     *  \param[in] subject : the subject (ie. the caller of the 'update' method).
     */
    virtual void update (EventInfo* evt, ISubject* subject) = 0;
};

/********************************************************************************/

/** \brief Interface of the Subject in the Observer Design Pattern.
 *
 *   This class provides methods for adding/removing observers to the subject.
 *
 *   There is also a 'notify' method that should be called when a notification has to be sent.
 *
 * Example:
 * \snippet observer1.cpp  snippet1
 *
 *  \see IObserver
 */
class ISubject
{
public:

    /** Destructor. */
    virtual ~ISubject() {}

    /** Returns the identifier the subject knows about.
     * \return identifier
     */
    virtual InterfaceId getInterface () = 0;

    /** Add an observer to the subject.
     * \param[in] observer : the observer to be added.
     */
    virtual void addObserver    (IObserver* observer) = 0;

    /** Remove an observer from the subject.
     * \param[in] observer : the observer to be removed.
     */
    virtual void removeObserver (IObserver* observer) = 0;

    /** Notify some piece of information to potential observer instances.
     * Actual implementations are likely to loop over all attached observer instances
     * and call the 'update' method on each of them.
     *
     * \param[in] event : information related to the notification
     */
    virtual void notify (EventInfo* event) = 0;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_DP_IOBSERVER_HPP_ */
