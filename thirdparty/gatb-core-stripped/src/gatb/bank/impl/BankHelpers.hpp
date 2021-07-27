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

/** \file BankHelpers.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Helpers for managing IBank objects
 */

#ifndef _GATB_CORE_BANK_IMPL_BANK_HELPERS_HPP_
#define _GATB_CORE_BANK_IMPL_BANK_HELPERS_HPP_

/********************************************************************************/

#include <gatb/bank/api/IBank.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/AbstractBank.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/api/IProperty.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
namespace impl      {
/********************************************************************************/

/** \brief Utility class holding useful methods for bank management
 */
class BankHelper
{
public:

    /** Singleton method.
     * \return the singleton instance. */
    static BankHelper& singleton()  { static BankHelper instance; return instance; }

    /** Convert one bank into another one.
     * \param[in] in  : the bank to be converted
     * \param[in] out : the converted bank
     * \param[in] progress : listener getting conversion progression information
     * */
    tools::misc::IProperties* convert (IBank& in, IBank& out, tools::dp::IteratorListener* progress=0);
};

/********************************************************************************/

/** \brief Bank implementation that delegates work to a referred bank.
 *
 * Implementation of the Proxy design pattern for the IBank interface.
 *
 * This class is not intended to be used by end users; it is rather used for
 * being subclassed.
 */
class BankDelegate : public AbstractBank
{
public:

    /** Constructor.
     * \param[in] ref : referred bank.
     */
    BankDelegate (IBank* ref) : _ref(0)  { setRef(ref); }

    /** Destructor. */
    ~BankDelegate () { setRef(0); }

    /** \copydoc AbstractBank::getId */
    std::string getId ()  { return _ref->getId(); }

	std::string getIdNb (int i)  { return _ref->getIdNb(i); }

	
    /** \copydoc AbstractBank::iterator */
    tools::dp::Iterator<Sequence>* iterator ()  { return _ref->iterator(); }

    /** \copydoc AbstractBank::getNbItems */
    int64_t getNbItems () { return _ref->getNbItems(); }

    /** \copydoc AbstractBank::insert */
    void insert (const Sequence& item)   { _ref->insert (item); }

    /** \copydoc AbstractBank::flush */
    void flush ()  {  _ref->flush ();  }

    /** \copydoc AbstractBank::getSize */
    u_int64_t getSize ()  { return _ref->getSize(); }

    /** \copydoc AbstractBank::estimate */
    void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize)  {  _ref->estimate (number, totalSize, maxSize);  }

    /** \copydoc AbstractBank::estimateNbItems */
    int64_t estimateNbItems () { return _ref->estimateNbItems(); }

	/** \copydoc AbstractBank::estimateNbItems */
	int64_t estimateNbItemsBanki (int i) { return _ref->estimateNbItemsBanki(i); }
	
    /** \copydoc AbstractBank::estimateSequencesSize */
    u_int64_t estimateSequencesSize ()  { return _ref->estimateSequencesSize(); }

    /** \copydoc AbstractBank::getEstimateThreshold */
    u_int64_t getEstimateThreshold ()  { return _ref->getEstimateThreshold(); }

    /** \copydoc AbstractBank::setEstimateThreshold */
    void setEstimateThreshold (u_int64_t nbSeq)  { _ref->setEstimateThreshold(nbSeq); }

protected:

    IBank* _ref;
    void setRef (IBank* ref)  { SP_SETATTR(ref); }
};

/********************************************************************************/

/** \brief Bank that can filter sequences through a provided functor.
 *
 * This is an utility class that allows to filter a referred bank with a functor.
 *
 * The functor must define the following method:
 * \code
 * bool operator() (const Sequence& seq)
 * \endcode
 *
 * -> true means that the sequence is iterated, false the sequence is filtered out.
 */
template<typename Filter> class BankFiltered : public BankDelegate
{
public:

    /** Constructor.
     * \param[in] ref : referred bank.
     * \param[in] filter : functor that filters sequence.
     */
    BankFiltered (IBank* ref, const Filter& filter) : BankDelegate (ref), _filter(filter)  {}

    /** \copydoc tools::collections::Iterable::iterator */
    tools::dp::Iterator<Sequence>* iterator ()
    {
        // We create one iterator from the reference
        tools::dp::Iterator<Sequence>* it = _ref->iterator ();

        // We get the composition for this iterator
        std::vector<tools::dp::Iterator<Sequence>*> iterators = it->getComposition();

        if (iterators.size() == 1)  { return new tools::dp::impl::FilterIterator<Sequence,Filter> (it, _filter); }
        else
        {
            // We are going to create a new CompositeIterator, we won't need the one we just got from the reference
            LOCAL(it);

            // We may have to encapsulate each sub iterator with the filter.
            for (size_t i=0; i<iterators.size(); i++)  {
            	iterators[i] = new tools::dp::impl::FilterIterator<Sequence,Filter> (iterators[i], _filter);
            }
            return new tools::dp::impl::CompositeIterator<Sequence> (iterators);
        }
    }

private:

    Filter _filter;
};

/********************************************************************************/

/* \brief Bank factory associated to the BankFiltered class
 */
template<typename Filter> class BankFilteredFactory : public IBankFactory
{
public:

    /** Constructor.
     * \param[in] delegateFormat : format of the delegate bank to be created
     * \param[in] filter : functor used to filtering out some sequences of the referred bank. */
    BankFilteredFactory (const std::string& delegateFormat, const Filter& filter) : _format(delegateFormat), _filter(filter)  {}

    /** \copydoc IBankFactory::createBank */
    IBank* createBank (const std::string& uri)
    {
        /** We create the reference bank. */
        IBank* ref = Bank::getFactory(_format)->createBank (uri);

        /** We encapsulate with a filtered bank. */
        return new BankFiltered<Filter> (ref, _filter);
    }

private:

    std::string _format;
    Filter      _filter;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IMPL_BANK_HELPERS_HPP_ */
