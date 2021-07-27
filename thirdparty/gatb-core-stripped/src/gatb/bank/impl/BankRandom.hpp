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

/** \file BankRandom.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Random bank format
 */

#ifndef _GATB_CORE_BANK_IMPL_BANK_RANDOM_HPP_
#define _GATB_CORE_BANK_IMPL_BANK_RANDOM_HPP_

/********************************************************************************/

#include <gatb/bank/impl/AbstractBank.hpp>

#include <vector>
#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
namespace impl      {
/********************************************************************************/

/** \brief Implementation of IBank for random banks
 *
 * This class generates random genomic data and can be used for test purpose.
 */
class BankRandom : public AbstractBank
{
public:

    /** Returns the name of the bank format. */
    static const char* name()  { return "random"; }

    /** Constructor.
     * \param[in] nbSequences : number of sequences of the random bank
     * \param[in] length : length of a sequence. */
    BankRandom (size_t nbSequences, size_t length);

    /** Destructor. */
    ~BankRandom ();

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
        Iterator(const BankRandom& bank);

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
        const BankRandom& _bank;
        int64_t   _rank;
        bool      _isDone;

        tools::misc::Data* _dataRef;
        void setDataRef (tools::misc::Data* dataRef)  { SP_SETATTR(dataRef); }
    };

protected:

    size_t _nbSequences;
    size_t _length;

    friend class Iterator;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IMPL_BANK_RANDOM_HPP_ */
