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

/** \file HostInfo.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief
 */

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_HOST_INFO_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_HOST_INFO_HPP_

/********************************************************************************/

#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/** \brief Host information
 */
class HostInfo
{
public:

    /** Get information about host
     * \return information as a IProperties instance
     */
    static IProperties& getInfo()
    {
        static system::SmartObject singleton;

        if (singleton.hasRef() == false)
        {
            IProperties* props = new Properties();

            props->add (0, "host");
            props->add (1, "name",             "%s",   system::impl::System::info().getHostName().c_str());
            props->add (1, "nb_cores",         "%d",   system::impl::System::info().getNbCores());
            props->add (1, "memory",           "%.1f", (double)system::impl::System::info().getMemoryPhysicalTotal() / (double)system::GBYTE);
            props->add (1, "disk_current_dir", "%.1f", (double)system::impl::System::file().getAvailableSpace(system::impl::System::file().getCurrentDirectory()) / (double)system::MBYTE);
            props->add (1, "max_file_nb",      "%lld", system::impl::System::file().getMaxFilesNumber());
            props->add (1, "pid",              "%d",   system::impl::System::thread().getProcess());

            singleton.setRef (props);
        }
        return * (dynamic_cast<IProperties*>(singleton.getRef()));
    }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_HOST_INFO_HPP_ */
