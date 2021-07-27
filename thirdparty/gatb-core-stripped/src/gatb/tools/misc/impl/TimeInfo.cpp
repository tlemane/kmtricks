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

#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/system/impl/System.hpp>

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
TimeInfo::TimeInfo () : _time(system::impl::System::time())
{
	_synchro = System::thread().newSynchronizer();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
TimeInfo::TimeInfo (system::ITime& aTime) : _time(aTime)
{
	_synchro = System::thread().newSynchronizer();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void TimeInfo::start (const char* name)
{
    _entriesT0 [name] = _time.getTimeStamp();
}

	
//destructor
TimeInfo::~TimeInfo()
{
	if(_synchro)
	{
		delete _synchro;
		_synchro = 0;
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
void TimeInfo::stop (const char* name)
{
    _entries [name] += _time.getTimeStamp() - _entriesT0 [name];
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
 ** REMARKS :
 GR: this must be thread safe ! because multiple threads may call this on the same TimeInfo object at the same time !
*********************************************************************/
TimeInfo& TimeInfo::operator+= (TimeInfo& ti)
{
	_synchro->lock();

	const std::map <std::string, u_int32_t>& entries = ti.getEntries();

    for (map <string, u_int32_t>::const_iterator it = entries.begin(); it != ti.getEntries().end(); ++it)
    {
        _entries[it->first] += it->second;
    }

	_synchro->unlock();
    return *this;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
 GR: this must be thread safe ! because multiple threads may call this on the same TimeInfo object at the same time !
*********************************************************************/
TimeInfo& TimeInfo::operator/= (size_t nb)
{
	_synchro->lock();

    for (map <string, u_int32_t>::const_iterator it = _entries.begin(); it != _entries.end(); ++it)
    {
        _entries[it->first] = (u_int32_t)  ((float)it->second / (float)nb);
    }
	
	_synchro->unlock();
    return *this;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
const std::map <std::string, u_int32_t>& TimeInfo::getEntries ()
{
    return _entries;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
u_int32_t TimeInfo::getEntryByKey (const std::string& key)
{
    u_int32_t result = 0;
    std::map <std::string, u_int32_t>::iterator  it = _entries.find (key);
    if (it != _entries.end())  { result = it->second; }
    return result;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
tools::misc::IProperties* TimeInfo::getProperties (const std::string& root)
{
    u_int32_t total = 0;

    /** We first compute the aggregated time. */
    std::map <std::string, u_int32_t>::const_iterator  it;
    for (it = getEntries().begin(); it != getEntries().end();  it++)   {  total += it->second;  }

    tools::misc::IProperties* props = new tools::misc::impl::Properties();

    props->add (0, root, "%.3f", (double)total/1000.0);

    for (it = getEntries().begin(); it != getEntries().end();  it++)
    {
        props->add (1, it->first.c_str(), "%.3f", (double)(it->second) / 1000.0);
    }

    return props;
}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/
