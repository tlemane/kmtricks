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

/** \file Cell.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Implementation of INode interface
 */

#ifndef _GATB_CORE_STORAGE_CELL_HPP_
#define _GATB_CORE_STORAGE_CELL_HPP_

#include <gatb/tools/storage/api/ICell.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace storage   {
namespace impl      {
/********************************************************************************/

/** \brief Partial implementation of the INode interface
 *
 * The 'remove' method is still abstract.
 */
class Cell : public virtual ICell, public system::SmartPointer
{
public:

    /** Constructor. */
    Cell (ICell* parent, const std::string& id)  : _parent(0), _id(id), _compressLevel(0)
	{
    	setParent(parent);
    	if (_parent != 0)  { _compressLevel = _parent->getCompressLevel(); }
	}

    /** Destructor. */
    ~Cell ()   {  setParent(0);  }

    /** \copydoc ICell::getParent  */
    ICell* getParent () const { return _parent; }

    /** \copydoc ICell::getId  */
    const std::string& getId ()  const { return _id; }

    /** \copydoc ICell::getFullId  */
    std::string getFullId (char sep='.') const
    {
        if (_parent != 0)   {  return _parent->getId().empty()==false ? _parent->getId() + sep + getId() : getId();  }
        else                {  return getId();  }
    }

    /** \copydoc ICell::setCompressLevel  */
    void setCompressLevel (int level)  { _compressLevel = level; }

    /** \copydoc ICell::getCompressLevel  */
    int getCompressLevel () const  { return _compressLevel; }

private:

    ICell* _parent;
    void setParent (ICell* parent)  {  _parent = parent;  }

    std::string _id;

    int _compressLevel;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_STORAGE_CELL_HPP_ */
