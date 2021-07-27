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

/** \file SystemInfoCommon.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Implementations common to various OS.
 */

#ifndef _GATB_CORE_SYSTEM_IMPL_SYSTEM_COMMON_HPP_
#define _GATB_CORE_SYSTEM_IMPL_SYSTEM_COMMON_HPP_

/********************************************************************************/

#include <gatb/system/api/ISystemInfo.hpp>
#include <gatb/system/api/IMemory.hpp>
#include <gatb/system/api/Exception.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace system    {
/** \brief Implementation of Operating System abstraction layer */
namespace impl      {
/********************************************************************************/

/** \brief Abstract class for ISystemInfo interface.
 *
 * This class factorizes some methods with implementations common to
 * several operating systems.
 *
 * It is still abstract since it doesn't implement all methods of
 * ISystemInfo
 */
class SystemInfoCommon : public ISystemInfo
{
public:

    /** \copydoc ISystemInfo::getVersion */
    std::string getVersion () const;

    /** \copydoc ISystemInfo::getBuildDate */
    std::string getBuildDate () const;

    /** \copydoc ISystemInfo::getBuildCompiler */
    std::string getBuildCompiler () const;

    /** \copydoc ISystemInfo::getBuildOptions */
    std::string getBuildOptions () const;

    /** \copydoc ISystemInfo::getBuildSystem */
    std::string getBuildSystem () const;

    /** \copydoc ISystemInfo::getHomeDirectory */
    std::string getHomeDirectory ()  const {  return getenv("HOME") ? getenv("HOME") : ".";  }
    
    /** \copydoc ISystemInfo::getMemoryPhysicalFree */
    u_int64_t getMemoryPhysicalFree () const  { return getMemoryPhysicalTotal()-getMemoryPhysicalUsed(); }

    /** \copydoc ISystemInfo::getMemoryProject */
    u_int64_t getMemoryProject () const  {  return std::min (getMemoryPhysicalFree() / (2*MBYTE), (u_int64_t)(5*1024)); }

    /** \copydoc ISystemInfo::getMemorySelfUsed */
    u_int64_t getMemorySelfUsed() const  { return 0; }

    /** \copydoc ISystemInfo::getMemorySelfUsed */
    u_int64_t getMemorySelfMaxUsed() const  { return 0; }

    /** \copydoc ISystemInfo::createCpuInfo */
    virtual CpuInfo* createCpuInfo (); //  { return new CpuInfoCommon(); }
};

/********************************************************************************/

/** \brief Linux implementation for ISystemInfo interface.
 */
class SystemInfoLinux : public SystemInfoCommon
{
public:

    /** \copydoc ISystemInfo::getNbCores */
    size_t getNbCores () const ;

    /** \copydoc ISystemInfo::getHostName */
    std::string getHostName () const ;

    /** \copydoc ISystemInfo::getMemoryPhysicalTotal */
    u_int64_t getMemoryPhysicalTotal () const ;

    /** \copydoc ISystemInfo::getMemoryPhysicalUsed */
    u_int64_t getMemoryPhysicalUsed () const ;

    /** \copydoc ISystemInfo::getMemoryBuffers */
    u_int64_t getMemoryBuffers () const ;

    /** \copydoc ISystemInfo::getMemorySelfUsed */
    u_int64_t getMemorySelfUsed() const;

    /** \copydoc ISystemInfo::getMemorySelfUsed */
    u_int64_t getMemorySelfMaxUsed() const;
};

/********************************************************************************/

/** \brief MacOs implementation for ISystemInfo interface.
 */
class SystemInfoMacos : public SystemInfoCommon
{
public:

    /** \copydoc ISystemInfo::getNbCores */
    size_t getNbCores () const ;

    /** \copydoc ISystemInfo::getHostName */
    std::string getHostName () const ;

    /** \copydoc ISystemInfo::getMemoryPhysicalTotal */
    u_int64_t getMemoryPhysicalTotal () const;

    /** \copydoc ISystemInfo::getMemoryPhysicalUsed */
    u_int64_t getMemoryPhysicalUsed () const;

    /** \copydoc ISystemInfo::getMemoryBuffers */
    u_int64_t getMemoryBuffers () const        { throw ExceptionNotImplemented(); }

    /** \copydoc ISystemInfo::getMemorySelfUsed */
    u_int64_t getMemorySelfUsed() const;
};

/********************************************************************************/

/** \brief Windows implementation for ISystemInfo interface.
 *
 * IMPORTANT ! Windows implementation is not fully done. Using it may lead to
 * ExceptionNotImplemented exception.
 */
class SystemInfoWindows : public SystemInfoCommon
{
public:

    /** \copydoc ISystemInfo::getNbCores */
    size_t getNbCores () const;

    /** \copydoc ISystemInfo::getHostName */
    std::string getHostName () const;

    /** \copydoc ISystemInfo::getMemoryPhysicalTotal */
    u_int64_t getMemoryPhysicalTotal () const  { throw ExceptionNotImplemented(); }

    /** \copydoc ISystemInfo::getMemoryPhysicalUsed */
    u_int64_t getMemoryPhysicalUsed ()  const   { throw ExceptionNotImplemented(); }

    /** \copydoc ISystemInfo::getMemoryBuffers */
    u_int64_t getMemoryBuffers ()       const  { throw ExceptionNotImplemented(); }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_SYSTEM_IMPL_SYSTEM_COMMON_HPP_ */
