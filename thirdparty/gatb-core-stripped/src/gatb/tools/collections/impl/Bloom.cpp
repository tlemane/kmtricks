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

#include <gatb/tools/collections/impl/Bloom.hpp>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/

u_int8_t bit_mask [] = {
    0x01,  //00000001
    0x02,  //00000010
    0x04,  //00000100
    0x08,  //00001000
    0x10,  //00010000
    0x20,  //00100000
    0x40,  //01000000
    0x80   //10000000
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

