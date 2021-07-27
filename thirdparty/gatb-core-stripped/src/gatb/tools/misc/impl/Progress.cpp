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

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/StringLine.hpp>
#include <gatb/system/impl/System.hpp>

#include <stdarg.h>
#include <stdio.h>

#define DEBUG(a)  //printf a

using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;

/********************************************************************************/
namespace gatb {  namespace core { namespace tools {  namespace misc {  namespace impl {
/********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Progress::Progress (u_int64_t ntasks, const char * msg, std::ostream& output)
    : os(output)
{
    message = StringLine::format (msg != NULL ? msg : "?");
    reset (ntasks);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Progress::reset (u_int64_t ntasks)
{
    todo    = ntasks;
    done    = 0;
    partial = 0;
    subdiv  = 100;
    steps   = (double)(ntasks>0 ? ntasks : 1) / (double)subdiv;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Progress::init ()
{
    /** We format the string to be displayed. */
    snprintf (buffer, 256, "| %-*s |\n",101, message.c_str());

    /** We dump the message.*/
    os << buffer;

    /** We forward the remaining of the init to postInit method. */
    postInit ();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Progress::finish ()
{
    /** We set the total done. */
    set (todo);

    /** We forward the remaining of the finish to postFinish method. */
    postFinish ();

    /** We reset all progression variables. */
    todo    = 0;
    done    = 0;
    partial = 0;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Progress::inc (u_int64_t ntasks_done)
{
    done    += ntasks_done;
    partial += ntasks_done;

    while (partial >= steps)
    {
        update (false);
        partial -= steps;
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
void Progress::set (u_int64_t ntasks_done)
{
    if (ntasks_done > done)
    {
        inc (ntasks_done-done);
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
void Progress::setMessage (const std::string& msg)
{
    message = StringLine::format (msg);

    update (false);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Progress::postInit ()
{
    os << "[";
    os.flush();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Progress::postFinish ()
{
    os << "]" << endl;
    os.flush();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Progress::update (bool first)
{
    os << "-";
    os.flush();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
ProgressTimer::ProgressTimer (u_int64_t ntasks, const char* msg, std::ostream& os)
    : Progress (ntasks, msg, os), heure_debut(0), heure_actuelle(0)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void ProgressTimer::postInit ()
{
    // os << endl;

    /** We get the current hour (in msec) */
    heure_debut = System::time().getTimeStamp();

    update (true);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void ProgressTimer::postFinish ()
{
    os << endl;
    os.flush();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void ProgressTimer::update (bool force)
{
    /** We get the current hour (in msec) */
    heure_actuelle = System::time().getTimeStamp();

    /** We compute the difference between the current time and the starting time. */
    double elapsed = (heure_actuelle - heure_debut) / 1000.0;

    /** A little check. */
    if (force || (done > 0))
    {
        fillBuffer (elapsed);

        /** We dump the string. */
        os << buffer;
        os.flush();
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
void ProgressTimer::fillBuffer (double elapsed)
{
    double speed  = elapsed > 0 ?  done / elapsed     : 0;
    double rem    = elapsed > 0 ? (todo-done) / speed : 0;

    if (done>todo) rem=0;

    int min_e  = (int)(elapsed / 60) ;
    elapsed -= min_e*60;
    int min_r  = (int)(rem / 60) ;
    rem -= min_r*60;

    /** We format the string to be displayed. */
    snprintf (buffer, sizeof(buffer), "%c[%s]  %-5.3g%%   elapsed: %3i min %-2.0f sec   remaining: %3i min %-2.0f sec",
        13,
        message.c_str(),
        100*(double)done/todo,
        min_e,elapsed,min_r,rem
    );
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
ProgressTimerAndSystem::ProgressTimerAndSystem (u_int64_t ntasks, const char* msg, std::ostream& os)
    : ProgressTimer (ntasks, msg, os), _cpuinfo(0), _memMax(0)
{
    setCpuInfo(System::info().createCpuInfo());
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
ProgressTimerAndSystem::~ProgressTimerAndSystem ()
{
    setCpuInfo (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void ProgressTimerAndSystem::fillBuffer (double elapsed)
{
    /** We get the memory used by the current process. */
    u_int64_t mem = System::info().getMemorySelfUsed() / 1024;
    if (_memMax<mem)  { _memMax=mem; }

    u_int64_t memMaxProcess = System::info().getMemorySelfMaxUsed() / 1024;

    /** We format the string to be displayed. */
    char tmp[128];
    snprintf (tmp, sizeof(tmp), "   cpu: %5.1f %%   mem: [%4lu, %4lu, %4lu] MB ",
        _cpuinfo->getUsage(),
        mem, _memMax, memMaxProcess
    );

    /** We call the parent method. */
    ProgressTimer::fillBuffer (elapsed);

    /** We concat to the final buffer. */
    strcat (buffer, tmp);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void ProgressTimerAndSystem::postInit ()
{
    _cpuinfo->start();
    _memMax = 0;

    ProgressTimer::postInit ();
}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/
