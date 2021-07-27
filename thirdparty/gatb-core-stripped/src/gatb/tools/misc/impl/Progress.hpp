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

/** \file Progress.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Progress feature
 */

#ifndef _GATB_CORE_TOOLS_MISC_PROGRESS_HPP_
#define _GATB_CORE_TOOLS_MISC_PROGRESS_HPP_

/********************************************************************************/

#include <gatb/system/api/types.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <string>
#include <iostream>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/** \brief Null implementation of IteratorListener.
 */
class ProgressNone : public dp::IteratorListener
{
public:
    /** Constructor. */
    ProgressNone (u_int64_t ntasks=0, const char* msg=0) {}
};

/********************************************************************************/

/** \brief Progress information display feature.
 *
 * This class is used to display some progression information, for instance during
 * the iteration of an iterator.
 *
 * It provides 3 methods that are supposed to be called by the progressing job:
 *      - init   : called at the beginning of the job
 *      - finish : called at the end of the job
 *      - inc    : called to notify some amount of job being done
 *
 * It can be used as this by the SubjectIterator class that encapsulates an iterator
 * and calls some functor for notifying progress information.
 *
 * This implementation merely dumps '-' characters as progression goes on.
 *
 * \see gatb::tools::dp::impl::SubjectIterator
 */
class Progress : public dp::IteratorListener
{
public:

    /** Constructor.
     * \param[in] ntasks : nb of items to be processed
     * \param[in] msg : message to be displayed
     * \param[in] os  : output stream where progress information is to be displayed
     */
    Progress (u_int64_t ntasks, const char * msg, std::ostream& os=std::cerr);

    /** Destructor. */
    virtual ~Progress () {}

    /** Initialization of the object. */
    void init ();

    /** Finish the progress information. */
    void finish ();

    /** Increase the number of currently done tasks.
     * \param[in] ntasks_done : amount of job done before previous call. */
    void inc (u_int64_t ntasks_done);

    /** Set the current number of tasks done.
     * \param[in] ntasks_done :  sets the current number of job done. */
    void set (u_int64_t ntasks_done);

    /** We reinitialize the progress object. */
    void reset (u_int64_t ntasks);

    /** \copydoc dp::IteratorListener::setMessage*/
    void setMessage (const std::string& msg);

protected:

    virtual void update     (bool force);
    virtual void postInit   ();
    virtual void postFinish ();

    std::string   message;
    u_int64_t     done;
    u_int64_t     todo;
    int           subdiv ; // progress printed every 1/subdiv of total to do
    double        partial;
    double        steps ; //steps = _todo/subidv
    std::ostream& os;
    char          buffer[512];

    friend class ProgressProxy;
};

/********************************************************************************/

/** \brief Progress information display feature with timing information
 *
 * This is another way of dumping information progressing, that provides information
 * such as the elapsed time and an estimation of the remaining time.
 */
class  ProgressTimer : public Progress
{
public:

    /** Constructor.
     * \param[in] ntasks : nb of items to be processed
     * \param[in] msg : message to be displayed
     * \param[in] os  : output stream where progress information is to be displayed
     */
    ProgressTimer (u_int64_t ntasks, const char* msg, std::ostream& os=std::cerr);

    /** Initialization of the object. */
    void init ()  { postInit (); }

protected:

    void update (bool force);
    void postInit ();
    void postFinish ();

    virtual void fillBuffer (double elapsed);

    system::ITime::Value  heure_debut;
    system::ITime::Value  heure_actuelle;
};

/********************************************************************************/

/** \brief Progress information display feature with timing and system information
 */
class  ProgressTimerAndSystem : public ProgressTimer
{
public:

    /** Constructor.
     * \param[in] ntasks : nb of items to be processed
     * \param[in] msg : message to be displayed
     * \param[in] os  : output stream where progress information is to be displayed
     */
    ProgressTimerAndSystem (u_int64_t ntasks, const char* msg, std::ostream& os=std::cerr);

    /** Destructor. */
    ~ProgressTimerAndSystem ();

protected:

    void postInit ();

    virtual void fillBuffer (double elapsed);

    system::ISystemInfo::CpuInfo* _cpuinfo;
    void setCpuInfo (system::ISystemInfo::CpuInfo* cpuinfo)  { SP_SETATTR(cpuinfo); }

    u_int64_t _memMax;
};

/********************************************************************************/

/* \brief Proxy for Progress class  */
class ProgressProxy : public dp::IteratorListener
{
public:

    ProgressProxy ()  : _ref(0) {}

    ProgressProxy (dp::IteratorListener* ref)  : _ref(0) { setRef(ref); }

    ProgressProxy (const ProgressProxy& p) : _ref(0) { setRef(p._ref); }

    ~ProgressProxy ()  { setRef(0); }

    /** Initialization of the object. */
    void init ()  { _ref->init(); }

    /** Finish the progress information. */
    void finish () { _ref->finish (); }

    /** Increase the number of currently done tasks.
     * \param[in] ntasks_done : amount of job done before previous call. */
    void inc (u_int64_t ntasks_done)  { _ref->inc (ntasks_done); }

    /** Set the current number of tasks done.
     * \param[in] ntasks_done :  sets the current number of job done. */
    void set (u_int64_t ntasks_done)  { _ref->set (ntasks_done); }

    /** */
    void reset (u_int64_t ntasks) { _ref->reset(ntasks); }

    /** \copydoc dp::IteratorListener::setMessage*/
    void setMessage (const std::string& msg)  { _ref->setMessage (msg); }

    dp::IteratorListener* getRef() const  { return _ref; }

private:
    dp::IteratorListener* _ref;
    void setRef (dp::IteratorListener* ref)  { SP_SETATTR(ref); }
};

/********************************************************************************/

/* \brief Synchro for Progress class  */
class ProgressSynchro : public ProgressProxy
{
public:

    ProgressSynchro () : _synchro(0)  {}

    ProgressSynchro (dp::IteratorListener* ref, system::ISynchronizer* synchro)
        : ProgressProxy (ref), _synchro(0)  { setSynchro(synchro); }

    ProgressSynchro (const ProgressSynchro& p) : ProgressProxy(p.getRef()), _synchro(0)  {  setSynchro(p._synchro); }

    ~ProgressSynchro ()  { setSynchro (0); }

    /** \copydoc dp::IteratorListener::init */
    void init ()  { system::LocalSynchronizer l(_synchro); ProgressProxy::init(); }

    /** \copydoc dp::IteratorListener::finish */
    void finish () { system::LocalSynchronizer l(_synchro); ProgressProxy::finish ();  }

    /** \copydoc dp::IteratorListener::inc*/
    void inc (u_int64_t ntasks_done)  { system::LocalSynchronizer l(_synchro);  ProgressProxy::inc (ntasks_done); }

    /** \copydoc dp::IteratorListener::set*/
    void set (u_int64_t ntasks_done)  { system::LocalSynchronizer l(_synchro);  ProgressProxy::set (ntasks_done); }

    /** */
    void reset (u_int64_t ntasks) { system::LocalSynchronizer l(_synchro);  ProgressProxy::reset(ntasks); }

    /** \copydoc dp::IteratorListener::setMessage*/
    void setMessage (const std::string& msg)  { system::LocalSynchronizer l(_synchro);  ProgressProxy::setMessage (msg); }

private:

    system::ISynchronizer* _synchro;
    void setSynchro (system::ISynchronizer* synchro)  { SP_SETATTR(synchro); }
};

/********************************************************************************/

/** We define a default class for progress information. */
typedef ProgressTimerAndSystem  ProgressDefault;

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_PROGRESS_HPP_ */
