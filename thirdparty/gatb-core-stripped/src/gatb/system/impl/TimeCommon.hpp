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

/** \file TimeCommon.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Implementations common to various OS.
 */

#ifndef _GATB_CORE_SYSTEM_IMPL_TIME_COMMON_HPP_
#define _GATB_CORE_SYSTEM_IMPL_TIME_COMMON_HPP_

/********************************************************************************/

#include <gatb/system/api/ITime.hpp>
#include <gatb/system/api/Exception.hpp>

#include <assert.h>
#include <sys/time.h>
#include <unistd.h>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace system    {
/** \brief Implementation of Operating System abstraction layer */
namespace impl      {
/********************************************************************************/

/** \brief Abstract implementation of ITime interface.
 *
 * This implementation provides a default getDateString method.
 */
class TimeAbstract : public ITime
{
public:

    /** Constructor.
     * \param[in] unit : unit used for this instance */
    TimeAbstract (Unit unit) : _unit (unit)  {}

    /** \copydoc ITime::getUnit */
    Unit getUnit ()   { return _unit; }

    /** \copydoc ITime::getUnit */
    std::string getDateString ()
    {
        char buffer[256];

        const char* format = "%4.4d%2.2d%2.2d_%2.2d%2.2d%2.2d";

        time_t now;       time (&now);
        struct tm today;  today = * localtime (&now);

        snprintf (buffer, sizeof(buffer), format,
            today.tm_year + 1900, today.tm_mon  + 1, today.tm_mday, today.tm_hour, today.tm_min, today.tm_sec
        );

        return buffer;
    }

protected:
    Unit _unit;
};

/********************************************************************************/

/** \brief Implementation of ITime interface.
 *
 * This class inherits from TimeAbstract class and add the getTimeStamp() method by using the
 * system gettimeofday function.
 */
class TimeSystem : public TimeAbstract
{
public:

    /** \copydoc TimeAbstract::TimeAbstract */
    TimeSystem (Unit unit) : TimeAbstract (unit)
    {
        /** We have to check that the provided unit is ok. */
        if (unit!=ITime::USEC  &&  unit!=ITime::MSEC  &&  unit!=ITime::SEC)
        {
            throw Exception ("bad time unit for TimeSystem");
        }
    }

    /** \copydoc ITime::getTimeStamp */
    Value getTimeStamp ()
    {
        timeval t;
        int err = gettimeofday (&t, NULL);

        if (err == 0)
        {
            switch (_unit)
            {
                case USEC:   return 1000000*t.tv_sec +  t.tv_usec;
                case MSEC:   return 1000*t.tv_sec +  t.tv_usec / 1000;
                case SEC:    return        t.tv_sec +  t.tv_usec / 1000000;
                default:     return 0;
            }
        }
        return 0;
    }
};

/********************************************************************************/

/** From Agner Fog's library. */
extern "C" u_int64_t ReadTSC ();

/**  \brief Implementation of ITime interface.
 *
 * This implementation tries to use a precise hardware time counter, known as
 * the Time Stamp Counter. It is then possible to get time stamps as number
 * of clock cycles which provides a great accuracy.
 *
 * Note that we use a third party library for this (see Agner Fog's work) providing
 * a assembler coded function that returns the time stamp counter value.
 */
class TimeCycle : public TimeAbstract
{
public:

    /** Constructor */
    TimeCycle () : TimeAbstract (UNDEFINED)  {}

    /** \copydoc ITime::getTimeStamp */
    Value getTimeStamp ()
    {
#if 0
        return ReadTSC();
#else
        return 0;
#endif
    }

    /** \return an estimation of the CPU clock frequency in GHz.*/
    double getClockFrequency ()
    {
        static double freq = 0;

        /** We compute the clock frequency at the first call. */
        if (freq == 0)
        {
            ITime::Value t0 = getTimeStamp();  sleep (1);  ITime::Value t1 = getTimeStamp();
            freq = (double) ((t1 - t0) / 1000000) / 1000.0;
        }
        return freq;
    }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_SYSTEM_IMPL_TIME_COMMON_HPP_ */
