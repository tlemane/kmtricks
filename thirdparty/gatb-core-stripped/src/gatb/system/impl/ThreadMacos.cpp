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

#ifdef __APPLE__

#include <gatb/system/impl/ThreadMacos.hpp>
#include <memory>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#include <unistd.h>

using namespace std;

/********************************************************************************/
namespace gatb { namespace core { namespace system { namespace impl {
/********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
class ThreadMacos : public IThread, public system::SmartPointer
{
public:
	ThreadMacos (void* (mainloop) (void*), void* data)  {
		
		
		//set stack size to 8 MB
		pthread_attr_t tattr;
		int ret = pthread_attr_init ( &tattr ) ;
		size_t size = 4096*2000 ; // must be multiple of page size
		ret = pthread_attr_setstacksize(&tattr, size);
		
		pthread_create (&_thread,&tattr,  mainloop, data);
		
		pthread_attr_destroy(&tattr);
		
	}
	
    ~ThreadMacos ()  { /* pthread_detach (_thread); */  }
    void join ()     { pthread_join   (_thread, NULL);  }
    Id getId ()  const { return (Id) _thread; }

private:
    pthread_t  _thread;
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
class SynchronizerMacos : public ISynchronizer, public system::SmartPointer
{
public:
	SynchronizerMacos ()            {  pthread_mutex_init (&_mutex, NULL);  }
    virtual ~SynchronizerMacos()    {  pthread_mutex_destroy (&_mutex);     }

    void   lock ()  { pthread_mutex_lock   (&_mutex); }
    void unlock ()  { pthread_mutex_unlock (&_mutex); }

private:
    pthread_mutex_t  _mutex;
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IThread* ThreadFactoryMacos::newThread (void* (*mainloop) (void*), void* data)
{
    return new ThreadMacos (mainloop, data);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
ISynchronizer* ThreadFactoryMacos::newSynchronizer (void)
{
    return new SynchronizerMacos ();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IThread::Id ThreadFactoryMacos::getThreadSelf()
{
    return (IThread::Id) pthread_self();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
u_int64_t ThreadFactoryMacos::getProcess ()
{
    return getpid ();
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* __LINUX__ */
