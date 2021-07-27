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

#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

using namespace std;
using namespace gatb::core::tools::dp;
using namespace gatb::core::system;

/********************************************************************************/
namespace gatb  {
namespace core  {
namespace tools {
namespace dp    {
namespace impl  {
/********************************************************************************/

/** */
class CommandStartSynchro : public ICommand, public system::SmartPointer
{
public:

    CommandStartSynchro (ICommand* ref, system::ISynchronizer* synchro) : _ref(0), _synchro(synchro)  { setRef(ref); }
    ~CommandStartSynchro ()  { setRef(0); }

    void execute ()
    {
        if (_ref && _synchro)
        {
            /** We lock/unlock the synchronizer. */
            _synchro->lock();  _synchro->unlock ();

            /** We execute the delegate command. */
            _ref->execute();
        }
    }

private:
    ICommand* _ref;
    void setRef (ICommand* ref)  { SP_SETATTR(ref); }

    system::ISynchronizer* _synchro;
};

/********************************************************************************/
class SynchronizerNull : public system::ISynchronizer, public system::SmartPointer
{
public:
    void   lock ()  {}
    void unlock ()  {}
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
size_t SerialDispatcher::dispatchCommands (std::vector<ICommand*>& commands, ICommand* post)
{
    TIME_START (ti, "compute");

    for (std::vector<ICommand*>::iterator it = commands.begin(); it != commands.end(); it++)
    {
        if (*it != 0)  {  (*it)->use ();  (*it)->execute ();  (*it)->forget ();  }
    }

    /** We may have to do some post treatment. Note that we do it in the current thread. */
    if (post != 0)  {  post->use ();  post->execute ();  post->forget ();  }

    TIME_STOP (ti, "compute");

    return ti.getEntryByKey("compute");
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
system::ISynchronizer* SerialDispatcher::newSynchro ()
{
    return new SynchronizerNull();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Dispatcher::Dispatcher (size_t nbUnits, size_t groupSize) : _nbUnits(nbUnits), _groupSize(groupSize)
{
    if (_nbUnits==0)  { _nbUnits = system::impl::System::info().getNbCores(); }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
size_t Dispatcher::dispatchCommands (std::vector<ICommand*>& commands, ICommand* postTreatment)
{
    TIME_START (ti, "compute");

    system::IThreadGroup* threadGroup = system::impl::ThreadGroup::create ();

    size_t idx = 0;

    /** We create threads and add them to the thread group. */
    for (std::vector<ICommand*>::iterator it = commands.begin(); it != commands.end(); it++, idx++)
    {
        /** We add the thread to the group. Note that we have to set two things:
         *  - the main loop function to be called by the thread
         *  - the information to be provided to the main loop.
         *
         *  Note that we provide here a IThreadGroup::Info instance as data for the main loop
         *  => it is important that at the beginning of the main loop we retrieved the correct
         *  type from the void* data. */
        threadGroup->add (
            mainloop,
            new IThreadGroup::Info (threadGroup,  new CommandStartSynchro (*it, threadGroup->getSynchro()), idx )
        );
    }

    /** We start the group. */
    threadGroup->start ();

    /** We may have to forward exceptions got in threads. */
    if (threadGroup->hasExceptions())  { throw threadGroup->getException(); }

    /** Some cleanup. */
    system::impl::ThreadGroup::destroy (threadGroup);

    TIME_STOP (ti, "compute");

    /** We return the result. */
    return ti.getEntryByKey("compute");
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
system::ISynchronizer* Dispatcher::newSynchro ()
{
    return system::impl::System::thread().newSynchronizer();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
system::IThread* Dispatcher::newThread (ICommand* command)
{
    return system::impl::System::thread().newThread (mainloop, command);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void* Dispatcher::mainloop (void* data)
{
    IThreadGroup::Info* info = (IThreadGroup::Info*) data;
    LOCAL (info);

    IThreadGroup* threadGroup = info->group;
    ICommand*    cmd          = (ICommand*) info->data;

    assert (threadGroup != 0);
    assert (cmd         != 0);

    /** Here, we should catch any exception thrown locally in a thread
     * and keep it to re-throw it in the main thread. */
    try
    {
        cmd->use ();
        cmd->execute();
        cmd->forget ();
    }
    catch (system::Exception& e)
    {
    	/** We copy the exception in the list of received exceptions of the thread group. */
        threadGroup->addException (e);
    }

    return 0;
}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/
