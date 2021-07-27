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

/** \file BankKmerIterator.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Kmer iterator on sequences of a bank
 */

#ifndef _GATB_CORE_KMER_IMPL_BANK_KMER_ITERATOR_HPP_
#define _GATB_CORE_KMER_IMPL_BANK_KMER_ITERATOR_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/bank/api/IBank.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

/* \brief Kmer iterator on data sequences from some bank.
 */
class BankKmerIterator :
    public tools::dp::Iterator<kmer_type>,
    public tools::dp::impl::AbstractSubjectIterator
{
public:

    /** Constructor.
     * \param[in] bank : the bank whose sequences are to be iterated.
     * \param[in] model : kmer model
     */
    BankKmerIterator (bank::IBank& bank, KmerModel& model)
        : _itSeq(0), _itKmer(model), _isDone(true),  _moduloMask(1), _current(0)
    {
        /** We create an iterator over the sequences of the provided bank.
         * Note that this is a dynamic allocation, so we will have to get rid of the instance
         * in the destructor. */
        setItSeq (bank.iterator());

        /** We set the modulo mask for which we will send notification to potential listeners.
         * Such notification is done on outer loop over Sequence objects. */
        _moduloMask = (1<<10) - 1;
    }

    /** Destructor. */
    ~BankKmerIterator ()
    {
        /** We get rid of the dynamically allocated Sequence iterator. */
        setItSeq(0);
    }

    /** \copydoc tools::dp::Iterator::first */
    void first()
    {
        /** We begin by notifying potential listeners that the iteration is beginning. */
        notifyInit ();

        /** We reset the iteration counter. We will check when this counter is equal to our modulo;
         * in such a case, we will notify our potential listeners.  */
        _current = 0;

        /** We go to the first item of the Sequence iteration. */
        _itSeq->first ();

        /** We use a shortcut variable in order to avoid some calls to the isDone method for the Sequence iteration.
         * This is important for performance concerns since the 'isDone' kmer iterator relies on it and that we use
         * a generic Iterator<Sequence>; in other words, we have here polymorphism on Sequence iterator and we have
         * to limit such polymorphic calls when the number of calls is huge (overhead due to polymorphism).
         */
        _isDone = _itSeq->isDone();

        /** We check that we have at least one sequence to iterate. */
        if (!_isDone)
        {
            /** We configure the kmer iterator with the current sequence data. */
            _itKmer.setData ((*_itSeq)->getData());

            /** We go to the first kmer. */
            _itKmer.first ();
        }
    }

    /** \copydoc tools::dp::Iterator::next */
    void next()
    {
        /** We look for the next kmer. */
        _itKmer.next ();

        /** We check the case where the kmer iteration is done. */
        if (_itKmer.isDone ())
        {
            /** We have no more kmer for the current sequence, therefore we go for the next sequence. */
            _itSeq->next();

            /** We check whether we have another sequence or not. */
            _isDone = _itSeq->isDone();
            if (!_isDone)
            {
                /** We configure the kmer iterator with the current sequence data. */
                _itKmer.setData ((*_itSeq)->getData());

                /** We go to the first kmer. */
                _itKmer.first ();

                /** We may have to notify potential listeners if we looped enough items. */
                if ((_current & _moduloMask) == 0)  { notifyInc (_current);  _current=0; }

                /** We increase the iterated kmers number. */
                _current++;
            }
        }
    }

    /** \copydoc tools::dp::Iterator::isDone */
    bool isDone()
    {
        /** If we are done, we notify potential listeners. */
        if (_isDone) { notifyFinish(); }

        /** We return the outer loop isDone status. */
        return _isDone;
    }

    /** \copydoc tools::dp::Iterator::item */
    kmer_type& item ()  { return _itKmer.item(); }

private:

    /** Outer loop iterator on Sequence. */
    tools::dp::Iterator<bank::Sequence>* _itSeq;
    void setItSeq (tools::dp::Iterator<bank::Sequence>* itSeq)  { SP_SETATTR(itSeq); }

    /** Inner loop iterator on kmer. */
    KmerModel::Iterator _itKmer;

    /** Shortcut (for performance). */
    bool _isDone;

    u_int32_t _current;
    u_int32_t _moduloMask;
};

/********************************************************************************/

/* \brief Kmer iterator on data sequences from some bank.
 */
template <typename kmer_type> class BankVectorKmerIterator :
        public tools::dp::Iterator<std::vector<kmer_type> >,
        public tools::dp::impl::AbstractSubjectIterator
{
public:

    /** Constructor.
     * \param[in] bank : the bank whose sequences are to be iterated.
     * \param[in] model : kmer model
     */
    BankVectorKmerIterator (bank::IBank& bank, KmerModel<kmer_type>& model)
        : _model(model), _itSeq(0), _isDone(true),  _moduloMask(1), _idx(0)
    {
        /** We create an iterator over the sequences of the provided bank.
         * Note that this is a dynamic allocation, so we will have to get rid of the instance
         * in the destructor. */
        setItSeq (bank.iterator());

        /** We set the modulo mask for which we will send notification to potential listeners.
         * Such notification is done on outer loop over Sequence objects. */
        _moduloMask = (1<<10) - 1;
    }

    /** Destructor. */
    ~BankVectorKmerIterator ()
    {
        /** We get rid of the dynamically allocated Sequence iterator. */
        setItSeq(0);
    }

    /** \copydoc tools::dp::Iterator::first */
    void first()
    {
        /** We begin by notifying potential listeners that the iteration is beginning. */
        notifyInit ();

        /** We reset the iteration counter. We will check when this counter is equal to our modulo;
         * in such a case, we will notify our potential listeners.  */
        _idx = 0;

        /** We go to the first item of the Sequence iteration. */
        _itSeq->first ();

        /** We use a shortcut variable in order to avoid some calls to the isDone method for the Sequence iteration.
         * This is important for performance concerns since the 'isDone' kmer iterator relies on it and that we use
         * a generic Iterator<Sequence>; in other words, we have here polymorphism on Sequence iterator and we have
         * to limit such polymorphic calls when the number of calls is huge (overhead due to polymorphism).
         */
        _isDone = _itSeq->isDone();

        /** We check that we have at least one sequence to iterate. */
        if (!_isDone)
        {
            _model.build ((*_itSeq)->getData(), *(this->_item));
        }
    }

    /** \copydoc tools::dp::Iterator::next */
    void next()
    {
		/** We have no more kmer for the current sequence, therefore we go for the next sequence. */
		_itSeq->next();

		/** We check whether we have another sequence or not. */
		_isDone = _itSeq->isDone();

		if (!_isDone)
		{
			/** We configure the kmer iterator with the current sequence data. */
		    _model.build ((*_itSeq)->getData(), *(this->_item));

			/** We may have to notify potential listeners if we looped enough items. */
			if ((_idx & _moduloMask) == 0)  { notifyInc (_idx);  _idx=0; }

			/** We increase the iterated kmers number. */
			_idx++;
		}
    }

    /** \copydoc tools::dp::Iterator::isDone */
    bool isDone()
    {
        /** If we are done, we notify potential listeners. */
        if (_isDone) { notifyFinish(); }

        /** We return the outer loop isDone status. */
        return _isDone;
    }

    /** \copydoc tools::dp::Iterator::item */
    std::vector<kmer_type>& item ()  { return *(this->_item); }

private:

    KmerModel<kmer_type>& _model;

    /** Outer loop iterator on Sequence. */
    tools::dp::Iterator<bank::Sequence>* _itSeq;
    void setItSeq (tools::dp::Iterator<bank::Sequence>* itSeq)  { SP_SETATTR(itSeq); }

    /** Shortcut (for performance). */
    bool _isDone;

    u_int32_t _idx;
    u_int32_t _moduloMask;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMPL_BANK_KMER_ITERATOR_HPP_ */
