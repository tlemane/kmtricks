/*****************************************************************************
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

/** \file BankZSTD.hpp
 *  \brief Zstandard-compressed FASTA/FASTQ bank format
 */

#ifndef _GATB_CORE_BANK_IMPL_BANK_ZSTD_HPP_
#define _GATB_CORE_BANK_IMPL_BANK_ZSTD_HPP_

/********************************************************************************/
#include <gatb/bank/impl/AbstractBank.hpp>

#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
namespace impl      {
/********************************************************************************/

/** \brief Implementation of IBank for zstd-compressed FASTA/FASTQ files
 *
 */
class BankZSTD : public AbstractBank
{
public:

    /** Returns the name of the bank format. */
    static const char* name()  { return "zstd"; }

    /** Constructor.
     * \param[in] filename : uri of the zstd-compressed file.
     */
    BankZSTD (const std::string& filename);

    /** Destructor. */
    ~BankZSTD ();

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

    /** \copydoc IBank::getSize
     */
    u_int64_t getSize ()  { return _filesize; }

    /** \copydoc IBank::estimate */
    void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize);

    /************************************************************/

    /** \brief Iterator for zstd-compressed FASTA/FASTQ files
     */
    class Iterator : public tools::dp::Iterator<Sequence>
    {
    public:

        /** Define what kind of comment we want to retrieve. Such comments can be
         * retrieved through gatb::core::bank::Sequence
         */
        enum CommentMode_e
        {
            /** Empty comments are provided to clients. */
            NONE,
            /** Comments with only the FASTA ID are provided to clients. */
            IDONLY,
            /** Full comments are provided to clients. */
            FULL
        };

        /** Constructor.
         * \param[in] ref : the associated iterable instance.
         * \param[in] commentMode : kind of comments we want to retrieve
         */
        Iterator (BankZSTD& ref, CommentMode_e commentMode = FULL);

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
        BankZSTD&    _ref;

        /** Tells what kind of comments we want as a client of the iterator. */
        CommentMode_e  _commentsMode;

        /** Tells whether the iteration is finished or not. */
        bool _isDone;

        /** Tells whether the instance is initialized. */
        bool _isInitialized;

        /** Number of times next has been called. */
        u_int64_t _nIters;

        /** Initialization method. */
        void init ();

        /** Finish method. */
        void finalize ();

        /** Decompressing stream (zstd_file_t defined in the .cpp). */
        void* _stream;

        /** Scratch strings used while parsing a record (buffered_strings_t). */
        void* _strings;

        bool get_next_seq (tools::misc::Vector<char>& data);
        bool get_next_seq (tools::misc::Vector<char>& data, std::string& comment, std::string& quality, CommentMode_e mode);

        size_t _index;
    };

protected:

    friend class Iterator;

    /** URI of the zstd file. */
    std::string _filename;

    /** Estimated uncompressed size. */
    u_int64_t _filesize;

    /** Initialization method (check existence, compute the file size). */
    void init ();
};

/********************************************************************************/

/** \brief Factory for the BankZSTD class. */
class BankZSTDFactory : public IBankFactory
{
public:

    /** \copydoc IBankFactory::createBank */
    IBank* createBank (const std::string& uri);
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IMPL_BANK_ZSTD_HPP_ */
