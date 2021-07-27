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

#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#define DEBUG(a)  printf a

using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

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
Algorithm::Algorithm (const std::string& name, int nbCores, gatb::core::tools::misc::IProperties* input)
    : _name(name), _input(0), _output(0), _info(0), _systemInfo(0), _dispatcher(0)
{
    setInput      (input ? input : new Properties());
    setOutput     (new Properties());
    setInfo       (new Properties());
    setSystemInfo (new Properties());

    if (nbCores < 0)  {  nbCores = _input->get(STR_NB_CORES)  ? _input->getInt(STR_NB_CORES) : 0;  }
    setDispatcher (new Dispatcher (nbCores) );

    _info->add (0, _name);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Algorithm::~Algorithm ()
{
    setInput      (0);
    setOutput     (0);
    setInfo       (0);
    setSystemInfo (0);
    setDispatcher (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Algorithm::run ()
{
    ISystemInfo::CpuInfo* cpuinfo = System::info().createCpuInfo();
    LOCAL (cpuinfo);

    cpuinfo->start();

    /** We execute the algorithm. */
    this->execute ();

    cpuinfo->stop();

    /** We gather some system information. */
    getSystemInfo()->add (1, "system");
    getSystemInfo()->add (2, "cpu",         "%.1f", cpuinfo->getUsage());
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
dp::IteratorListener* Algorithm::createIteratorListener (size_t nbIterations, const char* message)
{
    if (getInput()->get(STR_VERBOSE)==0)  { return new IteratorListener(); }

    switch (getInput()->getInt(STR_VERBOSE))
    {
        case 0: default:    return new IteratorListener ();
        case 1:             return new ProgressTimerAndSystem   (nbIterations, message);
        case 2:             return new ProgressTimer            (nbIterations, message);
        case 3:             return new Progress                 (nbIterations, message);
    }
}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/
