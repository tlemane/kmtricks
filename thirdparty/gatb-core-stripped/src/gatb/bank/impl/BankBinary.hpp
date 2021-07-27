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

/** \file BankBinary.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Binary bank format
 */

#ifndef _GATB_CORE_BANK__IMPL_BANK_BINARY_HPP_
#define _GATB_CORE_BANK__IMPL_BANK_BINARY_HPP_

/********************************************************************************/

#include <gatb/bank/impl/AbstractBank.hpp>

#include <vector>
#include <string>
#include <stdio.h>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
namespace impl      {
/********************************************************************************/

/** \brief Implementation of IBank for compressed format
 *
 * - a binary file is made of:
 *    - a magic number
 *    - a list of blocks
 *        - a block is:
 *              - one block size (on 4 bytes)
 *              - a list of sequences
 *                  - a sequence is:
 *                      - a sequence length (on 4 bytes)
 *                      - the nucleotides of the sequences (4 nucleotides encoded in 1 byte)
 * - number of sequences (on 4 bytes)
 *
 * Historically, BinaryBank has been used in the first step of the DSK tool to convert
 * one input FASTA file into a binary format. DSK used to read several times the reads
 * so having a binary (and so compressed) format had the nice effect to have less I/O
 * operations and therefore less execution time.
 *
 * In the following example, we can see how to convert any kind of bank into a binary bank:
 * \snippet bank8.cpp  snippet8_binary
 *
 */
class BankBinary : public AbstractBank
{
public:

    /** Returns the name of the bank format. */
    static const char* name()  { return "binary"; }

    /** Constructor. During a sequence insertion (see method 'insert'), a sequence may be split
     *  in sub sequences if invalid characters exist (like 'N'). A sub sequence is considered
     *  as valid if the number of consecutive letters is above some threshold (given as parameter).
     *  If this threshold is not provided, there is no split process during 'insert'
     * \param[in] filename : uri of the bank.
     * \param[in] nbValidLetters : threshold for sequence split in 'insert' method
     */
    BankBinary (const std::string& filename, size_t nbValidLetters=0);

    /** Destructor. */
    ~BankBinary ();

    /** \copydoc IBank::getId. */
    std::string getId ()  { return _filename; }

    /** \copydoc IBank::iterator */
    tools::dp::Iterator<Sequence>* iterator ()  { return new Iterator (*this); }

    /** \copydoc IBank::getNbItems */
    int64_t getNbItems () { return -1; }

    /** \copydoc IBank::insert */
    void insert (const Sequence& item);

    /** \copydoc IBank::flush */
    void flush ();

    /** \copydoc IBank::getSize */
    u_int64_t getSize ();

    /** \copydoc IBank::estimate */
    void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize);

    /** \copydoc IBank::remove. */
    void remove ();

    /** Set default buffer size (static method). 
      * \param[in] bufferSize : size of the buffer.    
      */
    static void setBufferSize (u_int64_t bufferSize);

    /** Check that the given uri is a correct binary bank. */
    static bool check (const std::string& uri);

    /************************************************************/

    /** \brief Specific Iterator impl for BankBinary class
     */
    class Iterator : public tools::dp::Iterator<Sequence>
    {
    public:
        /** Constructor.
         * \param[in] ref : the associated iterable instance.
         */
        Iterator (BankBinary& ref);

        /** Destructor */
        virtual ~Iterator ();

        /** \copydoc tools::dp::Iterator::first */
        void first();

        /** \copydoc tools::dp::Iterator::next */
        void next();

        /** \copydoc tools::dp::Iterator::isDone */
        bool isDone ()  { return _isDone; }

        /** \copydoc tools::dp::Iterator::item */
        Sequence& item ()
        {
            _item->getData().setEncoding (tools::misc::Data::BINARY);
            return *_item;
        }

        /** Estimation of the sequences information. */
        void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize);

    private:

        /** Reference to the underlying Iterable instance. */
        BankBinary&    _ref;

        /** Tells whether the iteration is finished or not. */
        bool _isDone;

        /** Block buffer read from file. */
        tools::misc::Data* _bufferData;
        void setBufferData (tools::misc::Data* bufferData)  { SP_SETATTR(bufferData); }

        int   cpt_buffer;
        int   blocksize_toread;
        int   nseq_lues;

        FILE* binary_read_file;

        size_t _index;
    };

protected:

    /** URI of the bank. */
    std::string _filename;

    size_t _nbValidLetters;

    unsigned char* buffer;
    int            cpt_buffer;
    int            read_write_buffer_size;
    FILE*          binary_read_file;

    void open  (bool write);
    void close ();
};

/********************************************************************************/

/* \brief Factory for the BankBinary class. */
class BankBinaryFactory : public IBankFactory
{
public:

    /** \copydoc IBankFactory::createBank */
    IBank* createBank (const std::string& uri)
    {
        return BankBinary::check(uri) ? new BankBinary (uri) : 0;
    }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK__IMPL_BANK_BINARY_HPP_ */
