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

/** \file ISystemInfo.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Interface providing information about the operating system
 */

#ifndef _GATB_CORE_SYSTEM_ISYSTEM_INFO_HPP_
#define _GATB_CORE_SYSTEM_ISYSTEM_INFO_HPP_

#include <gatb/system/api/types.hpp>
#include <gatb/system/api/ISmartPointer.hpp>
#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace system    {
/********************************************************************************/

/** \brief Interface providing some general information about the system.
 */
class ISystemInfo
{
public:

    /** Returns the version of the library.
     * \return the version. */
    virtual std::string getVersion () const = 0;

    /** Returns the date of the library generation.
     * \return the generation date. */
    virtual std::string getBuildDate () const = 0;

    /** Returns the compiler name
     * \return the compiler name. */
    virtual std::string getBuildCompiler () const = 0;

    /** Returns the compilation options
     * \return the compilation options. */
    virtual std::string getBuildOptions () const = 0;

    /** Returns the operating system name used for the library generation
     * \return the os name. */
    virtual std::string getBuildSystem () const = 0;

    /** Returns the number of available cores.
     * \return the number of cores. */
    virtual size_t getNbCores () const = 0;

    /** Returns the host name.
     * \return the host name. */
    virtual std::string getHostName () const = 0;

    /** Returns home directory.
     * \return the home directory uri. */
    virtual std::string getHomeDirectory () const = 0;

    /** Get the size (in bytes) of the physical memory
     * \return the physical memory size */
    virtual u_int64_t getMemoryPhysicalTotal () const = 0;

    /** Get the size (in bytes) of the used physical memory
     * \return the used physical memory size */
    virtual u_int64_t getMemoryPhysicalUsed () const = 0;

    /** Get the size (in bytes) of the free physical memory
     * \return the free physical memory size */
    virtual u_int64_t getMemoryPhysicalFree () const = 0;

    /** Get a memory size (NOTE: in MBytes) for executing a program.
     * It may be the whole physical memory, some part of it or a constant size.
     * \return the project memory size */
    virtual u_int64_t getMemoryProject () const = 0;

    /** Get the size (in bytes) of the buffers memory
     * \return the buffers memory size */
    virtual u_int64_t getMemoryBuffers () const = 0;

    /** Get the size (in KBytes) of the memory used by the current process
     * \return the memory value */
    virtual u_int64_t getMemorySelfUsed() const = 0;

    /** Get the size (in KBytes) of the max memory used by the current process
     * \return the memory value */
    virtual u_int64_t getMemorySelfMaxUsed() const = 0;

    /** Destructor. */
    virtual ~ISystemInfo ()  {}

    /********************************************************************************/

    /** \brief Interface providing a way to get CPU usage information
     */
    class CpuInfo : public SmartPointer
    {
    public:

        /** Start CPU information acquisition. */
        virtual void start () = 0;

        /** Stop CPU information acquisition. */
        virtual void stop () = 0;

        /** Get the CPU usage between start and stop. */
        virtual double getUsage() = 0;
    };

    /** Create a CpuInfo object
     * \return the created object. */
    virtual CpuInfo* createCpuInfo () = 0;
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_SYSTEM_ISYSTEM_INFO_HPP_ */
