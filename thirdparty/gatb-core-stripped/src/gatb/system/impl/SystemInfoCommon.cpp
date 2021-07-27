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

#include <gatb/system/impl/SystemInfoCommon.hpp>
#include <gatb/system/api/build_info.hpp>

#include <string.h>
#include <unistd.h>
#include <sys/times.h>
#include <iostream>

using namespace std;

/********************************************************************************/
namespace gatb { namespace core { namespace system { namespace impl {
/********************************************************************************/

/*********************************************************************
         #####   #######  #     #  #     #  #######  #     #
        #     #  #     #  ##   ##  ##   ##  #     #  ##    #
        #        #     #  # # # #  # # # #  #     #  # #   #
        #        #     #  #  #  #  #  #  #  #     #  #  #  #
        #        #     #  #     #  #     #  #     #  #   # #
        #     #  #     #  #     #  #     #  #     #  #    ##
         #####   #######  #     #  #     #  #######  #     #
*********************************************************************/

/********************************************************************************/
/** \brief Interface providing a way to get CPU usage information
 */
class CpuInfoCommon : public ISystemInfo::CpuInfo
{
public:

    /** Start CPU information acquisition. */
    virtual void start ()
    {
        struct tms timeSample;
        CPU0     = times (&timeSample);
        SysCPU0  = timeSample.tms_stime;
        UserCPU0 = timeSample.tms_utime;
    }

    /** Stop CPU information acquisition. */
    virtual void stop ()
    {
        struct tms timeSample;
        CPU1     = times (&timeSample);
        SysCPU1  = timeSample.tms_stime;
        UserCPU1 = timeSample.tms_utime;
    }

    /** Get the CPU usage between start and stop. */
    virtual double getUsage()
    {
        stop ();

        double percent = 0;

        if (CPU1 <= CPU0 || SysCPU1 < SysCPU0 ||  UserCPU1 < UserCPU0)
        {
            percent = -1.0;
        }
        else
        {
            percent = (SysCPU1 - SysCPU0) +  (UserCPU1 - UserCPU0);
            percent /= (CPU1 - CPU0);
            percent *= 100;
        }
        return percent;
    }

private:

    clock_t CPU0, SysCPU0, UserCPU0;
    clock_t CPU1, SysCPU1, UserCPU1;
};

/** */
ISystemInfo::CpuInfo* SystemInfoCommon::createCpuInfo ()  {  return new CpuInfoCommon (); }

std::string SystemInfoCommon::getVersion () const  { return STR_LIBRARY_VERSION; }

std::string SystemInfoCommon::getBuildDate () const { return STR_COMPILATION_DATE; }

std::string SystemInfoCommon::getBuildCompiler () const  { return STR_COMPILER; }

std::string SystemInfoCommon::getBuildOptions () const { return STR_COMPILATION_FLAGS; }

std::string SystemInfoCommon::getBuildSystem () const { return STR_OPERATING_SYSTEM; }

/*********************************************************************
                #        ###  #     #  #     #  #     #
                #         #   ##    #  #     #   #   #
                #         #   # #   #  #     #    # #
                #         #   #  #  #  #     #     #
                #         #   #   # #  #     #    # #
                #         #   #    ##  #     #   #   #
                #######  ###  #     #   #####   #     #
*********************************************************************/

#ifdef __linux__

#include "sys/sysinfo.h"
#include <unistd.h>
#include <sys/resource.h>
#include <sys/times.h>

/********************************************************************************/
size_t SystemInfoLinux::getNbCores () const
{
    size_t result = 0;

    /** We open the "/proc/cpuinfo" file. */
    FILE* file = fopen ("/proc/cpuinfo", "r");
    if (file)
    {
        char buffer[256];
        while (fgets(buffer, sizeof(buffer), file))  {  if (strstr(buffer, "processor") != NULL)  { result ++;  }  }
        fclose (file);
    }

    if (result==0)  { result = 1; }

    return result;
}

/********************************************************************************/
string SystemInfoLinux::getHostName () const
{
    string result;

    char hostname[1024];
    hostname[1023] = '\0';
    gethostname (hostname, sizeof(hostname)-1);
    result.assign (hostname, strlen(hostname));

    return result;
}

/********************************************************************************/
u_int64_t SystemInfoLinux::getMemoryPhysicalTotal () const
{
    struct sysinfo memInfo;
    sysinfo (&memInfo);
    u_int64_t result = memInfo.totalram;
    result *= memInfo.mem_unit;
    return result;
}

/********************************************************************************/
u_int64_t SystemInfoLinux::getMemoryPhysicalUsed () const
{
    struct sysinfo memInfo;
    sysinfo (&memInfo);
    u_int64_t result = memInfo.totalram - memInfo.freeram;
    result *= memInfo.mem_unit;
    return result;
}

/********************************************************************************/
u_int64_t SystemInfoLinux::getMemoryBuffers () const
{
    struct sysinfo memInfo;
    sysinfo (&memInfo);
    u_int64_t result = memInfo.bufferram;
    result *= memInfo.mem_unit;
    return result;
}

/********************************************************************************/
u_int64_t SystemInfoLinux::getMemorySelfUsed () const
{
    u_int64_t result = 0;
    FILE* file = fopen("/proc/self/status", "r");
    if (file)
    {
        char line[128];

        while (fgets(line, 128, file) != NULL)
        {
            if (strncmp(line, "VmRSS:", 6) == 0)
            {
                char* loop = line;
                result = strlen(line);
                while (*loop < '0' || *loop > '9') loop++;
                loop[result-3] = '\0';
                result = atoi(loop);
                break;
            }
        }
        fclose(file);
    }
    else
        std::cout << "warning: could not fopen /proc/self/status" << std::endl;
    return result;
}

/********************************************************************************/
u_int64_t SystemInfoLinux::getMemorySelfMaxUsed () const
{
    u_int64_t result = 0;
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage)==0)  {  result = usage.ru_maxrss;  }
    return result;
}
#endif

/*********************************************************************
            #     #     #      #####   #######   #####
            ##   ##    # #    #     #  #     #  #     #
            # # # #   #   #   #        #     #  #
            #  #  #  #     #  #        #     #   #####
            #     #  #######  #        #     #        #
            #     #  #     #  #     #  #     #  #     #
            #     #  #     #   #####   #######   #####
*********************************************************************/

#ifdef __APPLE__

#include <sys/types.h>
#include <sys/sysctl.h>

#include <mach/vm_statistics.h>
#include <mach/mach_types.h>
#include <mach/mach_init.h>
#include <mach/mach_host.h>
#include <mach/mach.h>

/********************************************************************************/
size_t SystemInfoMacos::getNbCores () const
{
    int numCPU = 0;

    int mib[4];
    size_t len = sizeof(numCPU);

    /* set the mib for hw.ncpu */
    mib[0] = CTL_HW;
    mib[1] = HW_AVAILCPU;  // alternatively, try HW_NCPU;

    /* get the number of CPUs from the system */
    sysctl (mib, 2, &numCPU, &len, NULL, 0);

    if (numCPU < 1)
    {
        mib[1] = HW_NCPU;
        sysctl (mib, 2, &numCPU, &len, NULL, 0 );
    }

    if (numCPU < 1)  {  numCPU = 1;  }

    return numCPU;
}

/********************************************************************************/
string SystemInfoMacos::getHostName () const
{
    string result;

    char hostname[1024];
    hostname[1023] = '\0';
    gethostname (hostname, sizeof(hostname)-1);
    result.assign (hostname, strlen(hostname));

    return result;
}

/********************************************************************************/
u_int64_t SystemInfoMacos::getMemoryPhysicalTotal () const
{
	int mib[2] = { CTL_HW, HW_MEMSIZE };
	size_t namelen = sizeof(mib) / sizeof(mib[0]);
	u_int64_t size;
	size_t len = sizeof(size);

	if (sysctl(mib, namelen, &size, &len, NULL, 0) < 0)
	{
		throw Exception ("Unable to get physical memory");
	}
	else
	{
		return size;
	}
}

/********************************************************************************/
u_int64_t SystemInfoMacos::getMemoryPhysicalUsed () const
{
	// see http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
	vm_size_t page_size;
	mach_port_t mach_port;
	mach_msg_type_number_t count;
	vm_statistics_data_t vm_stats;

	mach_port = mach_host_self();
	count = sizeof(vm_stats) / sizeof(natural_t);
	if (KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
	    KERN_SUCCESS == host_statistics(mach_port, HOST_VM_INFO,
	                                    (host_info_t)&vm_stats, &count))
	{
		int64_t myFreeMemory = (int64_t)vm_stats.free_count * (int64_t)page_size;

		int64_t used_memory = ((int64_t)vm_stats.active_count +
	                   (int64_t)vm_stats.inactive_count +
	                   (int64_t)vm_stats.wire_count) *  (int64_t)page_size;
		return myFreeMemory + used_memory;
	}
	else
	{
		throw Exception ("Unable to get free memory");
	}
}

/********************************************************************************/
u_int64_t SystemInfoMacos::getMemorySelfUsed () const
{
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

    if (KERN_SUCCESS != task_info (mach_task_self(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count))
    {
        return -1;
    }

    return t_info.resident_size / 1024;
}

#endif

/*********************************************************************
        #     #  ###  #     #  ######   #######  #     #   #####
        #  #  #   #   ##    #  #     #  #     #  #  #  #  #     #
        #  #  #   #   # #   #  #     #  #     #  #  #  #  #
        #  #  #   #   #  #  #  #     #  #     #  #  #  #   #####
        #  #  #   #   #   # #  #     #  #     #  #  #  #        #
        #  #  #   #   #    ##  #     #  #     #  #  #  #  #     #
         ## ##   ###  #     #  ######   #######   ## ##    #####
*********************************************************************/

#ifdef __WINDOWS__

#include <windows.h>

/********************************************************************************/
size_t SystemInfoWindows::getNbCores () const
{
    size_t result = 0;

    SYSTEM_INFO sysinfo;
    GetSystemInfo (&sysinfo);
    result = sysinfo.dwNumberOfProcessors;

    if (result==0)  { result = 1; }

    return result;
}

/********************************************************************************/
string SystemInfoWindows::getHostName () const
{
    string result;

    TCHAR  infoBuf[1024];
    DWORD  bufCharCount = sizeof(infoBuf)/sizeof(infoBuf[0]);

    if (GetComputerName( infoBuf, &bufCharCount ) )
    {
        result.assign (infoBuf, strlen(infoBuf));
    }

    return result;
}

#endif

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
