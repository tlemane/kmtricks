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

/** \file ThreadMacos.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Thread management for Macos
 */

#ifndef _GATB_CORE_SYSTEM_IMPL_MACOS_THREAD_HPP_
#define _GATB_CORE_SYSTEM_IMPL_MACOS_THREAD_HPP_

#include <gatb/system/api/IThread.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace system    {
namespace impl      {
/********************************************************************************/

/** \brief Factory that creates IThread instances.
 *
 *  Thread creation needs merely the main loop function that will be called.
 */
class ThreadFactoryMacos : public IThreadFactory
{
public:
    /** \copydoc IThreadFactory::newThread */
    IThread* newThread (void* (*mainloop) (void*), void* data);

    /** \copydoc IThreadFactory::newSynchronizer */
    ISynchronizer* newSynchronizer (void);

    /** \copydoc IThreadFactory::getThreadSelf */
    IThread::Id getThreadSelf();

    /** \copydoc IThreadFactory::getProcess */
    u_int64_t getProcess ();
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_SYSTEM_IMPL_MACOS_THREAD_HPP_ */
