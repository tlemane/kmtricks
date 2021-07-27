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

/** \file Range.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Iterable interface
 */

#ifndef _GATB_CORE_TOOLS_MISC_RANGE_HPP_
#define _GATB_CORE_TOOLS_MISC_RANGE_HPP_

/********************************************************************************/

#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/tools/collections/api/Iterable.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Tools package */
namespace tools     {
/** \brief Misc interfaces */
namespace misc      {
/********************************************************************************/

/* \brief Definition of an interval (inspired by std::pair). It is possible to define
 * 'reversed' range, ie. with a beginning greater than the end.
 */
template <class T> class Range : public collections::Iterable<T>, public system::SmartPointer
{
public:

    /** Default constructor. */
    Range() : begin(T()), end(T()) {}

    /** Constructor taking the [begin,end] couple as arguments.
     * \param[in] x : beginning of the range.
     * \param[in] y : end of the range.
     */
    Range(const T& x, const T& y) : begin(x), end(y) {}

    /** Copy constructor. */
    template <class U>  Range (const Range<U> &p) : begin(p.begin), end(p.end) { }

    /** \return the begin bound of the range. */
    T getBegin() const { return begin; }

    /** \return the end bound of the range. */
    T getEnd  () const { return end; }

    /** Returns the length of the range. */
    T getLength ()  const  { return (end >= begin ? end - begin + 1 : begin - end + 1); }

    /** Tells whether the provided value is included into the 'this' instance.
     * \param[in] val : value to be tested
     * \return true if the provided value is inside the current range, false otherwise.
     */
    bool includes (const T& val) const   {   return (this->begin <= val)  &&  (val <= this->end);  }

    /** Equality operator.
     * \param[in] r : the range to be compared to.
     * \return true if the ranges are the same (same beginning, same end), false otherwise. */
    bool operator== (const Range& r) const
    {
        return begin==r.begin && end==r.end;
    }

    /** InEquality operator.
     * \param[in] r : the range to be compared to.
     * \return false if the ranges are the same (same beginning, same end), true otherwise. */
    bool operator!= (const Range& r) const
    {
        return begin!=r.begin || end!=r.end;
    }

    /** */
    dp::Iterator<T>* iterator ()  { return new Iterator(*this); }

    /** */
    int64_t getNbItems ()       { return getLength(); }

    /** */
    int64_t estimateNbItems ()  { return getLength(); }

    /* */
    class Iterator : public dp::Iterator<T>
    {
    public:

        Iterator (const T& x, const T& y) : _begin(x), _end(y), _value(0), _isDone(true) {}

        Iterator (Range& ref) : _begin(ref.getBegin()), _end(ref.getEnd()), _value(0), _isDone(true)  {}

        void first()
        {
            _value = _begin;
            _isDone = _value > _end;
            if (!_isDone)  { *this->_item = _value; }
        }

        void next()
        {
            _value ++;
            _isDone = _value > _end;
            if (!_isDone)  { *this->_item = _value; }

        }

        bool isDone() { return _isDone; }

        T& item ()     { return (*this->_item); }

    private:
        T      _begin;
        T      _end;
        T      _value;
        bool   _isDone;
    };

    friend class Iterator;

private:
    /** The template (integer) type. */
    typedef T type;

    /** Beginning of the range. */
    T begin;

    /** End of the range. */
    T end;
};

/********************************************************************************/

/** We define a type for a range of kmer counts. */
typedef tools::misc::Range<CountNumber>  CountRange;

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_TOOLS_MISC_RANGE_HPP_ */
