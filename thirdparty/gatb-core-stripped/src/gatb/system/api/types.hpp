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

/** \file types.hpp
 *  \brief types definition for GATB.
 *  \date 01/03/2013
 *  \author edrezen
 *
 *   We define here some types used throughout the code.
 *
 *   Important: we define typedefs such as int16_t or u_int64_t. It is a good idea to use such typedefs
 *   instead of direct 'unsigned long' or 'short' for instance, because the actual number of used bytes
 *   may depend on the operating system/architecture. Using u_int32_t for instance ensure that we get
 *   an unsigned integer on 4 bytes.
 *
 *   Note that we use the <sys/types.h> file on Linux and MacOs. Such file may not exist on Windows (on Mingw
 *   to be more precise), so we propose here a definition. This is not perfect and should be improved.
 */

/********************************************************************************/

#ifndef _GATB_CORE_SYSTEM_TYPES_HPP_
#define _GATB_CORE_SYSTEM_TYPES_HPP_

/********************************************************************************/

#include <sys/types.h>
#include <vector>

/********************************************************************************/

/** We define a type for counting kmer occurrence. */
typedef int32_t CountNumber;

/** We define a type for a vector holding kmer counts. */
typedef std::vector<CountNumber> CountVector;

/********************************************************************************/

#endif /* _GATB_CORE_SYSTEM_TYPES_HPP_ */
