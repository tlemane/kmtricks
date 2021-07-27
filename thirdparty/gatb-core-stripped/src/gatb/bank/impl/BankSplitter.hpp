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

/** \file BankSplitter.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Interface definition for genomic databases management
 */

#ifndef _GATB_CORE_BANK_IMPL_BANK_SPLITTER_HPP_
#define _GATB_CORE_BANK_IMPL_BANK_SPLITTER_HPP_

/********************************************************************************/

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

/** \brief Split a given bank in pieces.
 *
 * This class allows to iterate the sequences of a given bank by splitting them.
 */
class BankSplitter : public AbstractBank
{
public:

    /** Returns the name of the bank format. */
    static const char* name()  { return "splitter"; }

    /** Constructor.
     * \param[in] reference : reference bank that provides the sequences to be split
     * \param[in] readMeanSize : size of a read
     * \param[in] overlap : nucleotides number that overlap between two consecutive reads
     * \param[in] coverage : number of occurrences for the reads
     */
    BankSplitter (
        IBank* reference,
        size_t readMeanSize,
        size_t   overlap,
        u_int8_t coverage
    );

    /** Destructor. */
    ~BankSplitter ();

    /** \copydoc IBank::getId. */
    std::string getId ()  { return "dummy"; }

    /** \copydoc IBank::iterator */
    tools::dp::Iterator<Sequence>* iterator ()  { return new Iterator (*this); }

    /** */
    int64_t getNbItems () { return -1; }

    /** \copydoc IBank::insert */
    void insert (const Sequence& item) {}

    /** \copydoc IBank::flush */
    void flush ()  {}

    /** \copydoc IBank::getSize */
    u_int64_t getSize ()  { return 0; }

    /** \copydoc IBank::estimate */
    void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize);

    /** \return maximum number of files. */
    static const size_t getMaxNbFiles ()  { return 0; }

    /************************************************************/

    class Iterator : public tools::dp::Iterator<Sequence>
    {
    public:

        /** */
        Iterator(const BankSplitter& bank);

        /** */
        ~Iterator();

        /** \copydoc tools::dp::Iterator::first */
        void first();

        /** \copydoc tools::dp::Iterator::next */
        void next();

        /** \copydoc tools::dp::Iterator::isDone */
        bool isDone ()  { return _isDone; }

        /** \copydoc tools::dp::Iterator::item */
        Sequence& item ()     { return *_item; }

    private:
        tools::misc::Data* _dataRef;
        void setDataRef (tools::misc::Data* dataRef)  { SP_SETATTR(dataRef); }

        tools::dp::Iterator<Sequence>* _itRef;
        void setItRef (tools::dp::Iterator<Sequence>* itRef)  { SP_SETATTR(itRef); }

        size_t    _readMeanSize;
        int64_t   _rank;
        int64_t   _nbMax;
        size_t    _offsetMax;
        size_t    _overlap;
        bool      _isDone;

        typedef std::pair<size_t,size_t> Offset;
        std::vector<Offset> _offsets;
    };

protected:

    IBank* _reference;
    void setReference (IBank* reference) { SP_SETATTR(reference); }

    size_t      _readMeanSize;
    u_int8_t    _coverage;
    size_t      _overlap;

    friend class Iterator;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IMPL_BANK_SPLITTER_HPP_ */
