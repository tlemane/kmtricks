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

#include <gatb/system/impl/System.hpp>

/********************************************************************************/
namespace gatb { namespace core { namespace system { namespace impl {
/********************************************************************************/
	
std::list<ThreadGroup*> ThreadGroup::_groups;

static pthread_mutex_t groupsMutex;
static int mutex_inited = 0;

#define LOCK()      pthread_mutex_lock   (&groupsMutex);
#define UNLOCK()    pthread_mutex_unlock (&groupsMutex);

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
ThreadGroup::ThreadGroup()
:  _startSynchro(0)
{
	init_mutex_if_needed();
	
    _startSynchro = System::thread().newSynchronizer();
    if (_startSynchro)  { _startSynchro->lock(); }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
ThreadGroup::~ThreadGroup()
{
    if (_startSynchro)  { delete _startSynchro; }

    /** We delete each thread. */
    for (std::vector<system::IThread*>::iterator it = _threads.begin(); it != _threads.end(); it++)
    {
        delete (*it);
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void ThreadGroup::add (void* (*mainloop) (void*), void* data)
{
    IThread* thr = System::thread().newThread (mainloop, data);

    _threads.push_back (thr);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void ThreadGroup::start ()
{
    /** We unlock the synchronizer => all threads of the group begin at the same time. */
    if (_startSynchro)  { _startSynchro->unlock(); }

    /** We join each thread. */
    for (std::vector<system::IThread*>::iterator it = _threads.begin(); it != _threads.end(); it++)
    {
        (*it)->join ();
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void ThreadGroup::init_mutex_if_needed()
{
	int  already_inited = __sync_fetch_and_add(&mutex_inited, 1);
	if(! already_inited)
	{
		pthread_mutex_init(&groupsMutex, NULL);
	}
	else
	{
        __sync_fetch_and_add(&mutex_inited, -1);
	}
}
	
/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IThreadGroup* ThreadGroup::create ()
{
	init_mutex_if_needed();

    LOCK();

    ThreadGroup* tg = new ThreadGroup;

    _groups.push_back(tg);

    UNLOCK();

    return tg;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void ThreadGroup::destroy (IThreadGroup* thr)
{
	init_mutex_if_needed();
	
    LOCK();

    /** We look for the provided thread group. */
    for (std::list<ThreadGroup*>::iterator it = _groups.begin(); it != _groups.end(); it++)
    {
        if (*it == thr)
        {
            _groups.erase (it);
            delete thr;
            break;
        }
    }
	
    UNLOCK();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IThreadGroup* ThreadGroup::find (IThread::Id id)
{
	init_mutex_if_needed();

	LOCK();

    /** We look for the provided thread group. */
    for (std::list<ThreadGroup*>::iterator it = _groups.begin(); it != _groups.end(); it++)
    {
        ThreadGroup* group = *it;

        for (std::vector<IThread*>::iterator itThread = group->_threads.begin();  itThread != group->_threads.end(); itThread++)
        {
            if ((*itThread)->getId() == id)  {  UNLOCK();  return group;   }
        }
    }
	
	UNLOCK();

    return 0;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool ThreadGroup::findThreadInfo (IThread::Id id, std::pair<IThread*,size_t>& info)
{
	init_mutex_if_needed();

	LOCK();

    /** We look for the provided thread group. */
    for (std::list<ThreadGroup*>::iterator it = _groups.begin(); it != _groups.end(); it++)
    {
        ThreadGroup* group = *it;

    	size_t idx=0;
        for (std::vector<IThread*>::iterator itThread = group->_threads.begin();  itThread != group->_threads.end(); itThread++, idx++)
        {
            if ((*itThread)->getId() == id)
            {
            	info.first  = *itThread;
            	info.second = idx;

            	UNLOCK();
            	return true;
            }
        }
    }

	UNLOCK();

    return false;
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
