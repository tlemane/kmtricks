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

/** \file BankFasta.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief FASTA bank format
 */

#ifndef _GATB_CORE_BANK_IMPL_BANK_FASTA_HPP_
#define _GATB_CORE_BANK_IMPL_BANK_FASTA_HPP_

/********************************************************************************/
#include <zlib.h>

#include <gatb/bank/impl/AbstractBank.hpp>

#include <vector>
#include <string>

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
/** \brief Implementation for genomic databases management. */
namespace impl      {
/********************************************************************************/

/** \brief Implementation of IBank for FASTA format
 *
 *  This class provides FASTA management in GATB.
 *
 *  Actually, it provides FASTA and FASTQ formats, both in uncompressed and gzip formats.
 *
 *  In case of FASTQ files, the iterated Sequence objects will provide quality information.
 *
 * Sample of use (note however that it is better to use Bank::open for opening a bank):
 * \snippet bank1.cpp  snippet1_bank
 */
class BankFasta : public AbstractBank
{
public:

    /** Returns the name of the bank format. */
    static const char* name()  { return "fasta"; }

    /** Constructor.
     * \param[in] filename : uri of the bank.
     * \param[in] output_fastq : tells whether the file is in fastq or not.
     * \param[in] output_gz: tells whether the file is gzipped or not
     */
    BankFasta (const std::string& filename, bool output_fastq = false, bool output_gz = false);

    /** Destructor. */
    ~BankFasta ();

    /** \copydoc IBank::getId. */
    std::string getId ()  { return _filenames[0]; }

    /** \copydoc IBank::iterator */
    tools::dp::Iterator<Sequence>* iterator ()  { return new Iterator (*this); }

    /** \copydoc IBank::getNbItems */
    int64_t getNbItems () { return -1; }

    /** \copydoc IBank::insert */
    void insert (const Sequence& item);

    /** \copydoc IBank::flush */
    void flush ();

    /** \copydoc IBank::getSize */
    u_int64_t getSize ()  { return filesizes; }

    /** \copydoc IBank::estimate */
    void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize);

    static void setDataLineSize (size_t len) { _dataLineSize = len; }
    static size_t getDataLineSize ()  { return _dataLineSize; }

    /** \copydoc IBank::finalize */
    void finalize ();

    /************************************************************/

    /** \brief Specific Iterator impl for Bank class
     *
     * This implementation relies on the initial code from Minia. It wraps the
     * Iterator interface with the Minia code.
     *
     * Note that we made some effort not to put implementation code
     * here in the header; see in particular buffered_file and buffered_strings
     * attributes whose type is void* (and not the implementation type defined in
     * the cpp file).
     *
     * Note the small improvement compared to Minia: it is possible to create an
     * Iterator that provides (or not) sequence comments, according to the corresponding
     * parameter given to the Iterator constructor.
     *
     *  <b>IMPROVEMENTS</b>:
     *  - in case we have several banks to read, we could have at one time only one stream opened on the currently
     *  iterated file. The current implementation opens all streams, which may be avoided.
     */
    class Iterator : public tools::dp::Iterator<Sequence>
    {
    public:

        /** Define what kind of comment we want to retrieve. Such comments can be retrieved through
         * gatb::core::bank::Sequence
         */
        enum CommentMode_e
        {
            /** Empty comments are provided to clients. */
            NONE,
            /** Comments with only the FASTA ID are provided to clients. \n
             *  Ex: 'ENSTTRP00000009639pep:novel'*/
            IDONLY,
            /** Full comments are provided to clients. \n
             *  Ex: 'ENSTTRP00000001236pep:novel genescaffold:turTru1:GeneScaffold_3311:182575:189152:1'*/
            FULL
        };

        /** Constructor.
         * \param[in] ref : the associated iterable instance.
         * \param[in] commentMode : kind of comments we want to retrieve
         */
        Iterator (BankFasta& ref, CommentMode_e commentMode = FULL);

        /** Destructor */
        ~Iterator ();

        /** \copydoc tools::dp::Iterator::first */
        void first();

        /** \copydoc tools::dp::Iterator::next */
        void next();

        /** \copydoc tools::dp::Iterator::isDone */
        bool isDone ()  { return _isDone; }

        /** \copydoc tools::dp::Iterator::item */
        Sequence& item ()     { return *_item; }

        /** Estimation of the sequences information */
        void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize);

    private:

        /** Reference to the underlying Iterable instance. */
        BankFasta&    _ref;

        /** Tells what kind of comments we want as a client of the iterator. */
        CommentMode_e  _commentsMode;

        /** Tells whether the iteration is finished or not. */
        bool _isDone;

        /** Tells whether the instance is initialized. */
        bool _isInitialized;

        /* Number of time next has been called   */
        u_int64_t   _nIters;
        
        /** Initialization method. */
        void init ();

        /** Finish method. */
        void finalize ();

        int   index_file; // index of current file

        void** buffered_file;
        void*  buffered_strings;   // variable_string_t *read, *dummy, *header;

        bool get_next_seq           (tools::misc::Vector<char>& data);
        bool get_next_seq           (tools::misc::Vector<char>& data, std::string& comment, std::string& quality, CommentMode_e mode);

        bool get_next_seq_from_file (tools::misc::Vector<char>& data, int file_id);
        bool get_next_seq_from_file (tools::misc::Vector<char>& data, std::string& comment, std::string& quality, int file_id, CommentMode_e mode);

        size_t _index;
    };

protected:

    /** \return maximum number of files. */
    static const size_t getMaxNbFiles ()  { return 1; }

    friend class Iterator;

    bool _output_fastq;
    bool _output_gz;
    
    /** List of URI of the banks. */
    std::vector<std::string> _filenames;

    u_int64_t filesizes;  // estimate of total size for all files
    size_t    nb_files;   // total nb of files

    /** File handle for inserting sequences into the bank. */
    FILE* _insertHandle;

    gzFile _gz_insertHandle;
    
    static size_t _dataLineSize;

    /** Initialization method (compute the file sizes). */
    void init ();
};

/********************************************************************************/

/* \brief Factory for the BankFasta class. */
class BankFastaFactory : public IBankFactory
{
public:

    /** \copydoc IBankFactory::createBank */
    IBank* createBank (const std::string& uri);
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IMPL_BANK_FASTA_HPP_ */
