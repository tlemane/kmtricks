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

/** \file Command.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Implementation of ICommand and IDispatcher
 */

#ifndef _GATB_CORE_DP_ITERATOR_IMPL_COMMAND_HPP_
#define _GATB_CORE_DP_ITERATOR_IMPL_COMMAND_HPP_

#include <gatb/system/impl/System.hpp>
#include <gatb/tools/designpattern/api/ICommand.hpp>
#include <gatb/tools/misc/api/Macros.hpp>

#include <list>

/********************************************************************************/
namespace gatb  {
namespace core  {
namespace tools {
namespace dp    {
namespace impl  {
/********************************************************************************/

/** \brief Launches commands in current thread
 *
 * A dispatcher that uses the calling thread, so no parallelization.
 *
 * This implementation can be useful to process ICommand instances in an serial way
 * when it is required, while keeping an uniform API (ie. call dispatchCommands)
 * for running ICommand instances.
 */
class SerialDispatcher : public IDispatcher
{
public:

    /** Destructor (defined because of presence of virtual methods). */
    virtual ~SerialDispatcher() {}

    /** \copydoc IDispatcher::dispatchCommands */
    size_t dispatchCommands (std::vector<ICommand*>& commands, ICommand* postTreatment=0);

    /** \copydoc IDispatcher::getExecutionUnitsNumber */
    size_t getExecutionUnitsNumber () { return 1; }

    /** \copydoc IDispatcher::setGroupSize */
    void setGroupSize (size_t groupSize)  { }

    /** \copydoc IDispatcher::getGroupSize */
    size_t getGroupSize () const  { return 1; }

private:

    /** */
    system::ISynchronizer* newSynchro ();
};

/********************************************************************************/

/** \brief Launches commands in different threads
 *
 *  This implementation launches commands through different threads => parallelization.
 *
 *  This implementation of IDispatcher is central in a tool design because it allows
 *  to uses all available cores. If one knows the number N of available cores on the computer,
 *  one has just to split some job by creating N ICommand instances and then just dispatch these
 *  commands through a Dispatcher: each command will be launched in a separated
 *  thread, and, thanks to the operating system architecture, each thread should be processed
 *  on an available core.
 *
 *  Note: it wouldn't be reasonable to use more ICommand instances than available cores.
 *  By default, if the number of dispatching units is not provided in the constructor of
 *  Dispatcher, it retrieves the number of available cores through the
 *  a call to system functions, and uses it as default value. This means
 *  that default constructor will use by default the whole CPU multicore power.
 */
class Dispatcher : public IDispatcher
{
public:

    /** Constructor.
     * \param[in] nbUnits : number of threads to be used. If 0 is provided, one tries to guess the number of available cores.
     * \param[in] groupSize : number of items to be retrieved from the iterator by one thread in a synchronized way
     */
    Dispatcher (size_t nbUnits=0, size_t groupSize=0);

    /** \copydoc IDispatcher::dispatchCommands */
    size_t dispatchCommands (std::vector<ICommand*>& commands, ICommand* postTreatment=0);

    /** \copydoc IDispatcher::getExecutionUnitsNumber */
    size_t getExecutionUnitsNumber () { return _nbUnits; }

    /** \copydoc IDispatcher::setGroupSize */
    void setGroupSize (size_t groupSize)  { _groupSize = groupSize; }

    /** \copydoc IDispatcher::getGroupSize */
    size_t getGroupSize () const  { return _groupSize; }

private:

    /** */
    system::ISynchronizer* newSynchro ();

    /** */
    system::IThread* newThread (ICommand* command);

    /** */
    static void* mainloop (void* data);

    /** Number of execution units to be used for command dispatching. */
    size_t _nbUnits;

    /** Group size */
    size_t _groupSize;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DP_ITERATOR_IMPL_COMMAND_HPP_ */
