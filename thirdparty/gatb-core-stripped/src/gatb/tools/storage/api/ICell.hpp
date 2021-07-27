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

/** \file ICell.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Node
 */

#ifndef _GATB_CORE_STORAGE_ICELL_HPP_
#define _GATB_CORE_STORAGE_ICELL_HPP_

/********************************************************************************/

#include <gatb/system/api/ISmartPointer.hpp>
#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace storage   {
/********************************************************************************/

/** \brief Interface to be used for storage feature.
 *
 * This interface defines a few methods for hierarchical management of entities
 * that can be stored (likely in file system).
 */
class ICell : public virtual system::ISmartPointer
{
public:

    /** Destructor. */
    virtual ~ICell() {}

    /** Get the parent node (if any).
     * \return the parent node. */
    virtual ICell* getParent () const = 0;

    /** Return the identifier of the node.
     * \return the id of the node. */
    virtual const std::string& getId () const = 0;

    /** Return the full identifier (like a path "x.y.z")
     * \param[in] sep : separator character for the path string
     * \return the full identifier. */
    virtual std::string getFullId (char sep='.') const = 0;

    /** Physically remove the node. */
    virtual void remove () = 0;

    /** Get the root of the given cell
     * \param[in] cell : the cell we want to get the root
     * \return the root of the given cell  */
    static ICell* getRoot (ICell* cell)
    {
        ICell* loop=0;  for (loop=cell ; loop && loop->getParent() != 0;  loop=loop->getParent())  {}
        return loop;
    }

    /** Set the compression level (if supported)
     * \param[in] level : from 0 (no compression) to 9 (best compression). */
    virtual void setCompressLevel (int level) = 0;

    /** Get the compression level
     * \return the compression level. */
    virtual int getCompressLevel () const = 0;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_STORAGE_ICELL_HPP_ */
