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

/** \file Exception.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Operating System common abstraction.
 */

#ifndef _GATB_CORE_SYSTEM_IRESOURCE_HPP_
#define _GATB_CORE_SYSTEM_IRESOURCE_HPP_

/********************************************************************************/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <string>
#include <sstream>
#include <list>

/********************************************************************************/
namespace gatb      {
/** \brief Core package of the GATP project.
 *
 * The gatb::core package holds all the fundamental packages needed for writing
 * assembly algorithms.
 *
 * It holds some generic tools, like operating system abstraction, collections management or design patterns
 * concerns. It also holds recurrent needs as reading genomic banks, handling kmers and so on.
 */
namespace core      {
/** \brief Operating System abstraction layer */
namespace system    {
/********************************************************************************/

 /** \brief Exception class for operating system failures
  *
  * This class is used in GATB-CORE for throwing typed exceptions.
  */
 class Exception
 {
 public:

     /** Default constructor. */
     Exception ()  {}

     /** Constructor with some information in a "printf" way.
      * \param[in] format : format of the message
      * \param[in] ... : arguments for the format parameter
      */
     Exception (const char* format, ...)
     {
         va_list args;  va_start (args, format);  init (format, args);  va_end (args);
     }

     /** Returns the description message.
      * \return the message
      */
     const char* getMessage () const  { return _message.c_str(); }

 protected:

     /** */
     void init (const char* format, va_list args)
     {
         char buffer[256];
         vsnprintf (buffer, sizeof(buffer), format, args);
         _message.assign (buffer);
     }

     /** The informative message. */
     std::string _message;
 };

 /********************************************************************************/

 /** \brief Composite exception
  *
  * This class allows to build one exception message from the messages of several
  * exceptions.
  *
  * It may be used for instance by the ThreadGroup class.
  */
 class ExceptionComposite : public Exception
 {
 public:

	 /** Constructor
	  * \param[in] exceptions : list of exceptions from which the composite message is built. */
     ExceptionComposite (const std::list<Exception>& exceptions)
     {
         std::stringstream ss;
         for (std::list<Exception>::const_iterator it = exceptions.begin(); it != exceptions.end(); it++)
         {
             ss << it->getMessage() << std::endl;
         }
         _message = ss.str();
     }
 };

 /********************************************************************************/

 /** \brief Exception class with information got from strerror_r
  */
 class ExceptionErrno : public Exception
 {
 public:

    /** Constructor. The error message is built by calling strerror function.
     * \param[in] format : printf-like prototype
     */
    ExceptionErrno (const char* format, ...)
    {
        va_list args;  va_start (args, format);  init (format, args);  va_end (args);

        char* buffer = (char*) malloc (BUFSIZ);
        if (buffer != NULL)
        {
            *buffer = 0;

#ifdef __CYGWIN__
            // strerror_r doesnt seem to be declared in cygwin
            // "The strerror_r() function is similar to strerror(), but is thread safe."
            strerror (errno, buffer, BUFSIZ);
#else

#if !defined(__linux__) || ((_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600) && !defined(_GNU_SOURCE)) || !defined(__GLIBC__) // XSI-Compliant strerror_r
            int ret_code = strerror_r (errno, buffer, BUFSIZ);
            const char* ret_buffer = buffer;
            if (ret_code == 0)
#else // GNU's strerror_r might return a static string instead of filling buffer
                const char* ret_buffer = strerror_r (errno, buffer, BUFSIZ);
            if (ret_buffer != NULL)
#endif
#endif

            {  _message += std::string(" (") + std::string(ret_buffer) + std::string(")");  }
            free(buffer);
        }
    }
 };

 /********************************************************************************/

 /** \brief Exception class for lack of implementation
  */
 class ExceptionNotImplemented : public Exception
 {
 public:

     /** Constructor. */
     ExceptionNotImplemented ()  {  _message = "NOT IMPLEMENTED";  }
 };

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_SYSTEM_IRESOURCE_HPP_ */
