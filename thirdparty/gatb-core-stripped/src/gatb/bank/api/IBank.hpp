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

/** \file IBank.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Interface definition for genomic databases management
 */

#ifndef _GATB_CORE_BANK_IBANK_HPP_
#define _GATB_CORE_BANK_IBANK_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/collections/api/Bag.hpp>
#include <gatb/bank/api/Sequence.hpp>

/********************************************************************************/
namespace gatb      {
/** \brief Core package of the GATP project.
 *
 * The gatb::core package holds all the fundamental packages needed for writting
 * assembly algorithms.
 *
 * It holds some generic tools, like operating system abstraction, collections management or design patterns
 * concerns. It also holds recurrent needs as reading genomic banks, handling kmers and so on.
 */
namespace core      {
/** \brief Package for genomic databases management. */
namespace bank      {
/********************************************************************************/

/** \brief Interface for what we need to read genomic databases.
 *
 * The IBank interface provides means to:
 *  - read sequences from some container.
 *  - insert sequences into some container.
 *
 * One typical implementation of this interface is a FASTA bank parser.
 *
 * The key point here is that clients should use this interface instead of specific
 * implementations, wherever they need to get nucleotides sequences from some container.
 *
 * By doing this, such clients will support many bank formats (according to the fact that
 * a IBank implementation provides such a service). They will also support "fake" banks, for
 * instance random generated sequences or any kind of synthetic data.
 *
 * From a design point of view, IBank is an Iterable (something we can get iterators from) and
 * a Bag (something we can insert items into); in both cases, the item type is Sequence
 *
 * Some kind of factory method exists for using IBank instances: this is the gatb::core::bank::impl::Bank
 * class. This class allows to get a IBank instance by providing an URI. The actual implementation
 * class of the IBank interface is computed by analyzing the URI string itself, or the content of the
 * file defined by the URI. So, tools developers should in general get IBank instances this way, so their
 * tool would support different kind of input bank formats.
 *
 * \see Sequence
 * \see IBankFactory
 * \see impl::Bank
 */
class IBank : public tools::collections::Iterable<Sequence>, public tools::collections::Bag<Sequence>
{
public:

    /** Get an unique identifier for the bank (could be the URI of a FASTA file for instance).
     * \return the identifier */
    virtual std::string getId () = 0;

	/** In case of a composite bank, return the id of bank i
	 * \return id of sub bank i */
	virtual std::string getIdNb (int i) = 0;

	
	/** In case of a composite bank,
	 * \return estimation of the number of sequences of sub bank i */
	virtual int64_t estimateNbItemsBanki (int i) = 0;
	
	/** Return the vector of  sub IBank objects (in case of bank composite), or a vector containing only the bank itself
	 * \return the IBank objects. */
	virtual const std::vector<IBank*> getBanks() const  = 0;
	
    /** \copydoc tools::collections::Iterable::iterator */
    virtual tools::dp::Iterator<Sequence>* iterator () = 0;

    /** \copydoc tools::collections::Bag::insert */
    virtual void insert (const Sequence& item) = 0;

    /** Return the size of the bank (comments + data)
     *
     * The returned value may be an approximation in some case. For instance, if we use
     * a zipped bank, an implementation may be not able to give accurate answer to the
     * size of the original file.
     *
     * \return the bank size in bytes.*/
    virtual u_int64_t getSize () = 0;

    /** In case of a composite bank, return the number of sub banks.
     * \return number of sub banks. */
    virtual size_t getCompositionNb () = 0;

    /** Give an estimation of sequences information in the bank.
     * \param[out] number : sequences number
     * \param[out] totalSize : sequences size (in bytes)
     * \param[out] maxSize : max size size (in bytes)
     */
    virtual void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize) = 0;

    /** Shortcut to 'estimate' method.
     * \return estimation of the number of sequences */
    virtual int64_t estimateNbItems () = 0;

    /** Shortcut to 'estimate' method.
     * \return estimation of the size of sequences */
    virtual u_int64_t estimateSequencesSize () = 0;

    /** \return the number of sequences read from the bank for computing estimated information */
    virtual u_int64_t getEstimateThreshold () = 0;

    /** Set the number of sequences read from the bank for computing estimated information
     * \param[in] nbSeq : the number of sequences to be read.*/
    virtual void setEstimateThreshold (u_int64_t nbSeq) = 0;

    /** Remove physically the bank. This method will have non-empty implementation for banks using
     * file system for instance. */
    virtual void remove () = 0;

    /** Method that may be called when the bank is done. It is called by BankFasta destructor for instance. It will close fclose() or something equivalent. You don't need to call this function yourself. */
    virtual void finalize () = 0;
};

/********************************************************************************/

/** \brief Factory for IBank.
 *
 * This interface provides a factory method that builds a IBank* instance given some
 * identifier.
 *
 * Such an identifier can be an uri (FASTA banks for instance), or any mechanism allowing
 * to retrieve enough information for creating instances of a specific IBank implementation.
 *
 * Actually, the gatb::core::bank::impl::Bank class relies on a list of registered IBankFactory
 * instances.
 */
class IBankFactory : public system::SmartPointer
{
public:

    /** Create an instance of IBank for a given uri.
     * \param[in] uri : the uri used for create the bank
     * \return the IBank instance. */
    virtual IBank* createBank (const std::string& uri) = 0;
};

/********************************************************************************/

/** We define a type for bank identifiers (as integers) when dealing with several banks.*/
typedef u_int16_t BankIdType;

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IBANK_HPP_ */
