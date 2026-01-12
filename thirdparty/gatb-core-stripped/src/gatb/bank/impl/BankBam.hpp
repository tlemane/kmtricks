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

/** \file BankBam.hpp
 *  \date 01/07/2026
 *  \author langhorst
 *  \brief BAM bank format
 */

#ifndef _GATB_CORE_BANK_IMPL_BANK_BAM_HPP_
#define _GATB_CORE_BANK_IMPL_BANK_BAM_HPP_

/********************************************************************************/
#include <zlib.h>

#include <gatb/bank/impl/AbstractBank.hpp>

#include <vector>
#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
namespace impl      {
/********************************************************************************/

/** \brief Implementation of IBank for BAM format
 *
 *  This class provides BAM (Binary Alignment/Map) format support in GATB.
 *
 *  BAM files are compressed binary files containing sequence reads (aligned or unaligned).
 *  This implementation extracts only the read name and sequence data, ignoring:
 *  - Quality scores (left empty, similar to FASTA)
 *  - Alignment information (mapping position, CIGAR, flags, etc.)
 *
 *  BAM files use BGZF (Blocked GNU Zip Format) compression, which is compatible
 *  with standard gzip but allows random access through block boundaries.
 *
 * Sample of use:
 * \code
 * IBank* bank = Bank::open("reads.bam");
 * Iterator<Sequence>* it = bank->iterator();
 * for (it->first(); !it->isDone(); it->next()) {
 *     Sequence& seq = it->item();
 *     // Process sequence...
 * }
 * delete it;
 * delete bank;
 * \endcode
 */
class BankBam : public AbstractBank
{
public:

    /** Returns the name of the bank format. */
    static const char* name()  { return "bam"; }

    /** Constructor.
     * \param[in] filename : uri of the BAM file.
     */
    BankBam (const std::string& filename);

    /** Destructor. */
    ~BankBam ();

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
    u_int64_t getSize ()  { return _filesize; }

    /** \copydoc IBank::estimate */
    void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize);

    /************************************************************/

    /** \brief Iterator for BAM files
     *
     * This iterator reads BAM files using BGZF decompression.
     * It extracts read name and sequence from each alignment record,
     * while ignoring quality scores and alignment-specific data.
     */
    class Iterator : public tools::dp::Iterator<Sequence>
    {
    public:

        /** Constructor.
         * \param[in] ref : the associated iterable instance.
         */
        Iterator (BankBam& ref);

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
        BankBam&    _ref;

        /** Tells whether the iteration is finished or not. */
        bool _isDone;

        /** Tells whether the instance is initialized. */
        bool _isInitialized;

        /** Number of iterations */
        u_int64_t   _nIters;
        
        /** Initialization method. */
        void init ();

        /** Finish method. */
        void finalize ();

        /** BGZF file handle */
        void* _bgzf;

        /** Read next BAM record and populate sequence data */
        bool get_next_seq (tools::misc::Vector<char>& data, std::string& comment);

        /** Current sequence index */
        size_t _index;
    };

protected:

    friend class Iterator;

    /** URI of the BAM file */
    std::string _filename;

    /** Estimated file size */
    u_int64_t _filesize;

    /** Initialization method (compute the file size). */
    void init ();
};

/********************************************************************************/

/** \brief Factory for the BankBam class. */
class BankBamFactory : public IBankFactory
{
public:

    /** \copydoc IBankFactory::createBank */
    IBank* createBank (const std::string& uri);
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IMPL_BANK_BAM_HPP_ */
