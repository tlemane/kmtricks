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

/** \file System.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Entry point class for accessing operating system operations
 */

#ifndef _GATB_CORE_SYSTEM_IMPL_SYSTEM_HPP_
#define _GATB_CORE_SYSTEM_IMPL_SYSTEM_HPP_

/********************************************************************************/

#include <gatb/system/api/config.hpp>
#include <gatb/system/impl/MemoryCommon.hpp>
#include <gatb/system/impl/TimeCommon.hpp>
#include <gatb/system/impl/SystemInfoCommon.hpp>
#include <gatb/system/impl/FileSystemLinux.hpp>
#include <gatb/system/impl/FileSystemMacos.hpp>
#include <gatb/system/impl/ThreadLinux.hpp>
#include <gatb/system/impl/ThreadMacos.hpp>
#include <map>
#include <functional>
#include <list>
#include <vector>
#include <iostream>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace system    {
/** \brief Implementation of Operating System abstraction layer */
namespace impl      {

/********************************************************************************/

 /** \brief Entry point providing access to operating system resources.
  *
  *  The IResource class provides a unique entry point for accessing different kinds
  *  of operating system resources (threads, time, file, etc...).
  *
  *  Normally, the client should use such an instance for getting OS resources
  *  instead of directly call system defendant functions. This is important because it
  *  will ease the build of client tools for different OS/architecture; the only thing
  *  to do is to use a specific instance of IResource for matching the correct OS.
  */
class System
 {
 public:

    /********************************************************************************/
    /** Access for info methods. */
    static ISystemInfo&  info  ()
    {
#ifdef __linux__
        static SystemInfoLinux instance;  return instance;
#endif

#ifdef __APPLE__
        static SystemInfoMacos instance;  return instance;
#endif

#ifdef __WINDOWS__
        #warning "TO BE DONE..."
#endif
    }

    /********************************************************************************/
     /** Access for time methods. */
     static ITime&  time ()
     {
#ifdef __linux__
        static TimeSystem instance (ITime::MSEC);  return instance;
#endif

#ifdef __APPLE__
        static TimeSystem instance (ITime::MSEC);  return instance;
#endif

#ifdef __WINDOWS__
        #warning "TO BE DONE..."
#endif
     }

     /********************************************************************************/
     /** Access for file methods. */
     static IFileSystem&     file    ()
     {
#ifdef __linux__
        static FileSystemLinux instance;  return instance;
#endif

#ifdef __APPLE__
        static FileSystemMacos instance;  return instance;
#endif

#ifdef __WINDOWS__
        #warning "TO BE DONE..."
#endif
     }

     /********************************************************************************/
     /** Access for memory methods. */
     static IMemory&         memory  ()
     {
#ifdef __linux__
        static MemoryCommon instance (allocator(), MemoryOperationsCommon::singleton());  return instance;
#endif

#ifdef __APPLE__
        static MemoryCommon instance (allocator(), MemoryOperationsCommon::singleton());  return instance;
#endif

#ifdef __WINDOWS__
        #warning "TO BE DONE..."
#endif
     }

     /********************************************************************************/
     /** Access for thread methods. */
     static IThreadFactory&  thread  ()
     {
#ifdef __linux__
        static ThreadFactoryLinux instance;  return instance;
#endif

#ifdef __APPLE__
        static ThreadFactoryMacos instance;  return instance;
#endif

#ifdef __WINDOWS__
        #warning "TO BE DONE..."
#endif
     }

 private:

     static IMemoryAllocator& allocator()
     {
#ifdef GATB_USE_CUSTOM_ALLOCATOR
         static MemorySizeStore instance (MemoryAllocatorStdlib::singleton());  return instance;
#else
         return MemoryAllocatorStdlib::singleton();
#endif
     }
};

/********************************************************************************/

/** \brief Implementation of IThreadGroup
 *
 */
class ThreadGroup : public IThreadGroup, public system::SmartPointer
{
public:

	/** Create a IThreadGroup instance
	 * \return the IThreadGroup instance */
    static IThreadGroup* create ();
	
	/** Destroy a IThreadGroup */
    static void destroy (IThreadGroup* thr);

    /** Get a thread of the group from one id
     * \param[in] id : the thread id
     * \return the IThreadGroup instance if found, NULL otherwise */
    static IThreadGroup* find (IThread::Id id);

    /** Get thread information (IThread and index within the group);
     * \param[in] id : the thread id
     * \return true if found */
    static bool findThreadInfo (IThread::Id id, std::pair<IThread*,size_t>& info);

    /** \copydoc IThreadGroup::add */
    void add (void* (*mainloop) (void*), void* data);

    /** \copydoc IThreadGroup::start */
    void start ();

    /** \copydoc IThreadGroup::getSynchro */
    ISynchronizer* getSynchro()  { return _startSynchro; }

    /** \copydoc IThreadGroup::size */
    size_t size() const { return _threads.size(); }

    /** \copydoc IThreadGroup::operator[] */
    IThread* operator[] (size_t idx)  { return _threads[idx]; }

    /** \copydoc IThreadGroup::addException */
    void addException (system::Exception e)
    {
        LocalSynchronizer synchro (_startSynchro);
        _exceptions.push_back (e);
    }

    /** \copydoc IThreadGroup::hasExceptions */
    bool hasExceptions() const { return _exceptions.empty() == false; }

    /** \copydoc IThreadGroup::getException */
    Exception getException () const  { return ExceptionComposite(_exceptions); }

private:

    ThreadGroup ();
    ~ThreadGroup();

    /** Initialize a synchronizer at first call. */
	static void init_mutex_if_needed ();

    std::vector<IThread*>  _threads;
    system::ISynchronizer* _startSynchro;

    static std::list<ThreadGroup*> _groups;

    std::list<system::Exception> _exceptions;
};

/********************************************************************************/

/** \brief Facility to share a common resource between several threads.
 *
 * When using multithreading, one has to take care about reads/writes on a resource T
 * by several threads at the same time.
 *
 * There are two ways to cope with this:
 * 	- using some synchronization mechanism (see ISynchronizer); each time the resource T
 * 	is accessed, one has to lock the synchronizer before the access and unlock the
 * 	synchronizer after the access
 * 	- using local information in each thread instead of accessing the T shared resource;
 * 	at the end of the threads, one has to use all the gathered local information to update
 * 	the T resource in a serial way.
 *
 * 	The ThreadObject class provides a generic way to follow the second approach. One provides
 * 	the shared resource to it, then the resource is cloned in N threads and each cloned resource
 * 	can be locally accessed within its thread without needing synchronization process.
 * 	At the end of the threads, the ThreadObject provides a way to get each cloned local resource
 * 	and so the user can update the original resource with the gathered information.
 *
 * 	It is possible to give no initial resource to a ThreadObject instance; in such a case,
 * 	a default original resource is created by using the default constructor of the type of
 * 	the resource.
 *
 * 	Sample of use:
 * \snippet multithreading5.cpp  snippet5_threadobject
 *
 */
template<typename T> class ThreadObject
{
public:

	/** Constructor
	 * \param[in] object : object to be shared by several threads.
	 * 		If none, default constructor of type T is used. */
    ThreadObject (const T& object = T()) : _object(object), _isInit(false), _synchro(0)
    {
        _synchro = system::impl::System::thread().newSynchronizer();
    }

    /** Destructor. */
    ~ThreadObject()
    {
        if (_synchro)  { delete _synchro; }

        for (typename std::map<IThread::Id,T*>::iterator it = _map.begin(); it != _map.end(); it++)  { delete it->second; }
    }

    /** Get the local T object for the current calling thread.
     * \return the local T object. */
    T& operator () ()
    {
        if (_isInit == false)
        {
            LocalSynchronizer ls (_synchro);
            if (_isInit == false)
            {
                /** We look for the TreadGroup if any. */
                IThreadGroup* group = ThreadGroup::find (System::thread().getThreadSelf());

                if (group)
                {
                    for (size_t i=0; i<group->size(); i++)
                    {
                        IThread* thread = (*group)[i];
                        T* newObject = new T(this->_object);
                        this->_map[thread->getId()] = newObject;
                        this->_vec.push_back(newObject);
                    }
                }
                else
                {
                	throw Exception ("ThreadObject::operator() : no defined group");
                }
                _isInit = true;
            }
        }

        return *(_map[System::thread().getThreadSelf()]);
    }

    /** Iterate all the local objects through a functor.
     * \param[in] fct : functor called for each local object. The functor must take a T or a const T& argument
     */
    template<typename Functor> void foreach (const Functor& fct)
    {  for (typename std::map<IThread::Id,T*>::iterator it = _map.begin(); it != _map.end(); it++)  { fct (*(it->second)); } }

    /** Get a pointer on the shared object.
     * \return a pointer on the shared object. */
    T* operator-> ()  { return &_object; }

    /** Get a reference on the shared object.
     * \return a reference on the shared object. */
    T& operator* ()  { return _object; }

    /** Get the number of local T objects
     * \return the number of local objects */
    size_t size() const { return _vec.size(); }

    /** Get the ith local T object
     * \param[in] idx : index of the local object to be returned
     * \return a reference on the wanted local object*/
    T& operator[] (size_t idx) { return *(_vec[idx]); }

private:
    std::map<IThread::Id,T*> _map;
    T _object;

    std::vector<T*> _vec;

    bool _isInit;
    system::ISynchronizer* _synchro;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#ifdef GATB_USE_CUSTOM_ALLOCATOR
    #define MALLOC      gatb::core::system::impl::System::memory().malloc
    #define CALLOC      gatb::core::system::impl::System::memory().calloc
    #define REALLOC     gatb::core::system::impl::System::memory().realloc
    #define FREE        gatb::core::system::impl::System::memory().free

    inline void* operator new      (size_t s)   { return MALLOC(s); }
    inline void* operator new[]    (size_t s)   { return MALLOC(s); }
    inline void  operator delete   (void* ptr)  { FREE(ptr); }
    inline void  operator delete[] (void* ptr ) { FREE(ptr); }

#else
    #define MALLOC      malloc
    #define CALLOC      calloc
    #define REALLOC     realloc
    #define FREE        free
#endif

#endif /* _GATB_SYSTEM_IMPL_SYSTEM_HPP_ */
