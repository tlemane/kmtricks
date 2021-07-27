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

/** \file BankConverterAlgorithm.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Bank conversion from one IBank to another IBank
 */

#ifndef _BANK_CONVERTER_ALGORITHM_HPP_
#define _BANK_CONVERTER_ALGORITHM_HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/bank/api/IBank.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/misc/api/Enums.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
namespace impl      {
/********************************************************************************/

/** \brief Algorithm implementation that converts a bank into a binary one
 *
 * This algorithm converts an input IBank instance into a binary output bank.
 *
 * It subclasses gatb::core::tools::misc::impl::Algorithm in order to get all the features
 * of its parent class.
 *
 * Historically, it was used by DSK to work on a binary bank instead of a FASTA
 * bank (ie. less I/O operations because of smaller size of binary banks)
 */
class BankConverterAlgorithm : public gatb::core::tools::misc::impl::Algorithm
{
public:

    /** Constructor.
     * \param[in] bank : bank to be converted (likely in FASTA format)
     * \param[in] kmerSize : kmer size
     * \param[in] outputUri : uri of the output binary bank. */
    BankConverterAlgorithm (IBank* bank, size_t kmerSize, const std::string& outputUri);

    /** Constructor. Used only to retrieved statistics/information gathered during
     * a previous execution of a BankConverterAlgorithm instance
     * \param[in] storage : storage instance where we get information from. */
    BankConverterAlgorithm (tools::storage::impl::Storage& storage);

    /** Destructor. */
    ~BankConverterAlgorithm ();

    /** \copydoc gatb::core::tools::misc::impl::Algorithm::execute */
    void execute ();

    /** Return the output binary bank
     * \return the IBank instance */
    IBank* getResult ()  { return _bankOutput; }

private:

    tools::misc::BankConvertKind _kind;

    IBank*      _bankInput;
    void setBankInput (IBank* bankInput)  { SP_SETATTR(bankInput); }

    IBank*      _bankOutput;
    void setBankOutput (IBank* bankOutput)  { SP_SETATTR(bankOutput); }

    std::string _outputUri;

    bank::IBank* createBank (
        tools::dp::Iterator<bank::Sequence>* inputSequences,
        size_t nbInputSequences,
        const std::string& outputName,
        u_int64_t& nbSeq,
        u_int64_t& sizeSeq
    );

    size_t _kmerSize;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _BANK_CONVERTER_ALGORITHM_HPP_ */

