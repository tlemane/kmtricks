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

/** \file AbstractBank.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Abstract implementation of the IBank interface
 */

#ifndef _GATB_CORE_BANK_IMPL_ABSTRACT_BANK_HPP_
#define _GATB_CORE_BANK_IMPL_ABSTRACT_BANK_HPP_

/********************************************************************************/

#include <gatb/bank/api/IBank.hpp>
#include <vector>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
namespace impl      {
/********************************************************************************/

/** \brief Abstract implementation of IBank for factorizing common behavior.
 *
 * This abstract implementation of the IBank interface provides some methods having the
 * same behavior for most implementations.
 *
 * Note that it implements the system::ISmartPointer interface as well, so it can be
 * used as a smart pointer.
 */
class AbstractBank : public IBank, public system::SmartPointer
{
public:

    /** Constructor. */
    AbstractBank () : _estimateThreshold(50000) {}

	
	std::string getIdNb (int i)  { return std::string("not_a_compo_bank"); }

	
	int64_t estimateNbItemsBanki (int i)  { return this->estimateNbItems(); }

	/** \copydoc IBank::getBanks */
	const std::vector<IBank*> getBanks() const  {
		std::vector<IBank*> _banks;
		_banks.push_back((IBank *)this);
		return _banks;
	};

	
    /** \copydoc IBank::estimateNbItems */
    int64_t estimateNbItems ()
    {
        u_int64_t number, totalSize, maxSize;    estimate (number, totalSize, maxSize);  return number;
    }

    /** \copydoc IBank::estimateSequencesSize */
    u_int64_t estimateSequencesSize ()
    {
        u_int64_t number, totalSize, maxSize;    estimate (number, totalSize, maxSize);  return totalSize;
    }

    /** \copydoc IBank::getEstimateThreshold */
    u_int64_t getEstimateThreshold ()  { return _estimateThreshold; }

    /** \copydoc IBank::setEstimateThreshold */
    void setEstimateThreshold (u_int64_t nbSeq) { _estimateThreshold = nbSeq; }

    /** \copydoc IBank::remove */
    void remove () {}

    /** \copydoc IBank::finalize */
    void finalize ()  {}

    /** \copydoc IBank::getCompositionNb */
    size_t getCompositionNb ()
    {
        tools::dp::Iterator<Sequence>* it = this->iterator();  LOCAL(it);
        return it->getComposition().size();
    }

private:

    u_int64_t _estimateThreshold;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IMPL_ABSTRACT_BANK_HPP_ */
