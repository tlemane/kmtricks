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

/** \file BankComposite.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Composite bank, ie. a bank made of other banks
 */

#ifndef _GATB_CORE_BANK_IMPL_BANK_COMPOSITE_HPP_
#define _GATB_CORE_BANK_IMPL_BANK_COMPOSITE_HPP_

/********************************************************************************/

#include <gatb/bank/impl/AbstractBank.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <vector>
#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
namespace impl      {
/********************************************************************************/

/** \brief IBank implementation for composite banks
 *
 * This class implements the Composite design pattern and allows to handle compound
 * banks the same way as a single one.
 *
 * This is mainly a base class for other subclasses (see BankAlbum) and is not often
 * instantiated directly.
 *
 * Most of the methods of IBank are implemented by iterating the list of referred banks.
 *
 * A BankComposite can be created with a vector of IBank instances to be referred.
 */
class BankComposite : public AbstractBank
{
public:

    /** Returns the name of the bank format. */
    static const char* name()  { return "composite"; }

    /** Default constructor (no referred bank). */
    BankComposite () : _nbItems(0), _size(0)  {}

    /** Constructor.
     * \param[in] banks : list of banks to be associated to the album.
     */
    BankComposite (const std::vector<IBank*>& banks) : _nbItems(0), _size(0)
    {
        for (size_t i=0; i<banks.size(); i++)  {  this->addBank (banks[i]); }
    }

    /** Destructor. */
    virtual ~BankComposite ()
    {
        for (size_t i=0; i<_banks.size(); i++)  { _banks[i]->forget(); }
    }

    /** \copydoc IBank::getId. */
    std::string getId ()
    {
        if (_id.empty())  {  for (size_t i=0; i<_banks.size(); i++)  {  if (i>0) { _id += ","; }  _id += _banks[i]->getId();  }  }
        return _id;
    }

	
	std::string getIdNb (int i)
	{
		if(i<(int)_banks.size())
			return _banks[i]->getId();
		else
			return
				this->getId();
	}
	
	
	int64_t estimateNbItemsBanki (int i)  {
		if(i<(int)_banks.size())
			return _banks[i]->estimateNbItems();
		else
			return  this->estimateNbItems();
	
	}

	
    /** Add a bank into the composite
     * \param[in] bank : the bank to be added. */
    void addBank (IBank* bank)  { bank->use();  _banks.push_back(bank); }

    /** \copydoc IBank::iterator */
    tools::dp::Iterator<Sequence>* iterator ()
    {
        std::vector <tools::dp::Iterator<Sequence>*>  iterators;
        for (size_t i=0; i<_banks.size(); i++)  { iterators.push_back (_banks[i]->iterator()); }
        return new tools::dp::impl::CompositeIterator<Sequence> (iterators);
    }

    /** \copydoc IBank::getNbItems */
    int64_t getNbItems ()
    {
        if (_nbItems == 0)  {  for (size_t i=0; i<_banks.size(); i++)  { _nbItems += _banks[i]->getNbItems(); } }
        return _nbItems;
    }

    /** \copydoc IBank::insert */
    void insert (const Sequence& item)   {  throw system::Exception ("Can't insert sequence in a composite bank.");  }

    /** \copydoc IBank::flush */
    void flush ()   {  for (size_t i=0; i<_banks.size(); i++)  { _banks[i]->flush(); } }

    /** \copydoc IBank::getSize */
    u_int64_t getSize ()
    {
        if (_size==0)  {   for (size_t i=0; i<_banks.size(); i++)  { _size += _banks[i]->getSize(); } }
        return _size;
    }

    /** \copydoc IBank::estimate */
    void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize)
    {
        number = totalSize = maxSize = 0;

        u_int64_t numberIth=0, totalSizeIth=0, maxSizeIth=0;
        for (size_t i=0; i<_banks.size(); i++)
        {
            _banks[i]->estimate (numberIth, totalSizeIth, maxSizeIth);
            number += numberIth;  totalSize += totalSizeIth;  maxSize = std::max (maxSize, maxSizeIth);
        }
    }

    /** */
    size_t getCompositionNb() { return _banks.size(); }

    /** \return maximum number of files. */
    static const size_t getMaxNbFiles ()  { return 30; }

    /** Return the vector of IBank objects.
     * \return the IBank objects. */
    const std::vector<IBank*> getBanks() const { return _banks; }

    /** Get the number of referred banks.
     * \return the number of referred banks */
    size_t getNbBanks () const { return _banks.size();  }

    /** Direct iteration of the IBank instances.
     * \param[in] fct : functor that will be called with each referred IBank object. */
    template<typename Functor> void iterateBanks (Functor fct)  {  for (size_t i=0; i<_banks.size(); i++)  { fct (*_banks[i], i); }  }

    /** Returns an iterator on the IBank objects (heap allocation)
     * \return the IBank iterator. */
    tools::dp::Iterator<IBank*>* iteratorBanks ()  { return new tools::dp::impl::VectorIterator<IBank*> (_banks); }

    /** \copydoc IBank::remove. */
    void remove ()  {  for (size_t i=0; i<_banks.size(); i++)  { _banks[i]->remove(); } }

protected:

    /** List of the banks. */
    std::vector<IBank*> _banks;

    u_int64_t _nbItems;
    u_int64_t _size;
    std::string _id;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IMPL_BANK_COMPOSITE_HPP_ */
