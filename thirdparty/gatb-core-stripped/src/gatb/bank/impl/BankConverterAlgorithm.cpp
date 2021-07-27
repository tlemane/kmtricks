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

#include <gatb/bank/impl/BankConverterAlgorithm.hpp>

#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/Banks.hpp>

#include <gatb/kmer/impl/Model.hpp>

// We use the required packages
using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

#define DEBUG(a)  //printf a

/********************************************************************************/
namespace gatb {  namespace core { namespace bank { namespace impl {
/********************************************************************************/

static const char* progressFormat1 = "Bank: fasta to binary                  ";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankConverterAlgorithm::BankConverterAlgorithm (IBank* bank,  size_t kmerSize, const std::string& outputUri)
: Algorithm ("bankconverter"), _kind(BANK_CONVERT_TMP), _bankInput(0), _bankOutput(0), _outputUri(outputUri), _kmerSize(kmerSize)
{
    setBankInput  (bank);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankConverterAlgorithm::BankConverterAlgorithm (tools::storage::impl::Storage& storage)
: Algorithm ("bankconverter"), _kind(BANK_CONVERT_NONE), _bankInput(0), _bankOutput(0), _kmerSize(0)
{
    string xmlString = storage(this->getName()).getProperty ("xml");
    stringstream ss; ss << xmlString;   getInfo()->readXML (ss);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankConverterAlgorithm::~BankConverterAlgorithm ()
{
    setBankInput  (0);
    setBankOutput (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankConverterAlgorithm::execute ()
{
    /** We may have no conversion at all to do. */
    if (_kind == BANK_CONVERT_NONE)
    {
        setBankOutput (_bankInput);
        return;
    }

    /** We may have to delete the binary file if it already exists. */
    System::file().remove (_outputUri);

    /** We get information about the FASTA bank. */
    u_int64_t number, totalSize, maxSize;
    _bankInput->estimate (number, totalSize, maxSize);

    /** We create the sequence iterator. */
    Iterator<Sequence>* itSeq = _bankInput->iterator();
    LOCAL (itSeq);

    u_int64_t   nbSeq = 0;
    u_int64_t sizeSeq = 0;

    /** We get the composition of the provided sequence iterator. */
    std::vector<Iterator<Sequence>*> iters = itSeq->getComposition();

    /** We use a block for measuring the time elapsed in it. */
    {
        TIME_INFO (getTimeInfo(), "conversion");

        /** We get information about the FASTA bank. */
        u_int64_t number, totalSize, maxSize;
        _bankInput->estimate (number, totalSize, maxSize);

        if (iters.size() == 1)
        {
            /** We set the output bank. */
            setBankOutput (createBank (itSeq, number, _outputUri, nbSeq, sizeSeq));
        }

        else if (iters.size() > 1)
        {
            vector<IBank*> ouputBanks;

            for (size_t i=0; i<iters.size(); i++)
            {
                /** We set the name of the current binary bank. */
                stringstream ss;   ss << _outputUri << i;

                /** We create a new binary bank and add it to the vector of output banks. */
                ouputBanks.push_back (createBank (iters[i], number / iters.size(), ss.str(), nbSeq, sizeSeq));
            }

            /** We set the result output bank. */
            setBankOutput (new BankAlbum (_outputUri, ouputBanks));
        }
        else
        {
            throw Exception ("Error in BankConverterAlgorithm : nbIterators=%d", iters.size());
        }
    }

    /** We flush the output bank (important if it is a file output). */
    _bankOutput->flush ();

    /** We gather some statistics. */
    getInfo()->add (1, "info");
    getInfo()->add (2, "input",            "%s",   _bankInput->getId().c_str());
    getInfo()->add (2, "composite_number", "%d",   iters.size());
    getInfo()->add (2, "sequences_number", "%ld",  nbSeq);
    getInfo()->add (2, "sequences_size",   "%ld",  sizeSeq);
    getInfo()->add (2, "output_size",      "%ld",  _bankOutput->getSize());
    getInfo()->add (2, "ratio",            "%.3f",  (double)sizeSeq / (double)_bankOutput->getSize());
    getInfo()->add (1, getTimeInfo().getProperties("time"));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IBank* BankConverterAlgorithm::createBank (
    Iterator<Sequence>* inputSequences,
    size_t              nbInputSequences,
    const string& outputName,
    u_int64_t& nbSeq,
    u_int64_t& sizeSeq
)
{
    DEBUG (("BankConverterAlgorithm::createBank  nbInputSequences=%ld  outputName='%s'  nbSeq=%d \n",
        nbInputSequences, outputName.c_str(), nbSeq
    ));

    /** We create a new binary bank. */
    IBank* result = new BankBinary (outputName, _kmerSize);

    /** We need an iterator on the input bank. */
    Iterator<Sequence>* itBank = createIterator<Sequence> (
        inputSequences,
        nbInputSequences,
        progressFormat1
    );
    LOCAL (itBank);

    /** We iterate the sequences of the input bank. */
    for (itBank->first(); !itBank->isDone(); itBank->next())
    {
        nbSeq ++;
        sizeSeq += (*itBank)->getDataSize();

        /** We insert the current sequence into the output bank. */
        result->insert (itBank->item());
    }

    return result;
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
