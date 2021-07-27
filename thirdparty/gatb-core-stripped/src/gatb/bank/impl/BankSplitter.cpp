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

#include <gatb/bank/impl/BankSplitter.hpp>
#include <gatb/tools/misc/api/StringsRepository.hpp>

#include <gatb/system/impl/System.hpp>

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;
using namespace gatb::core::tools::misc;

#define DEBUG(a)  //printf a

/********************************************************************************/
namespace gatb {  namespace core {  namespace bank {  namespace impl {
/********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankSplitter::BankSplitter (
    IBank*   reference,
    size_t   readMeanSize,
    size_t   overlap,
    u_int8_t coverage
)
    : _reference (0), _readMeanSize(readMeanSize), _coverage(coverage), _overlap(overlap)
{
    setReference (reference);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankSplitter::~BankSplitter ()
{
    setReference (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankSplitter::estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize)
{
    gatb::core::tools::dp::Iterator<gatb::core::bank::Sequence>* itSeq = _reference->iterator();
    LOCAL (itSeq);

    itSeq->first();
    assert (itSeq->isDone() == false);

    Data& data = itSeq->item().getData();

    size_t offsetMax = data.size() - _readMeanSize;
    size_t delta     = _readMeanSize -_overlap;
    size_t nb        = 1 + offsetMax / delta;

    number    = nb * _coverage;
    totalSize = number * _readMeanSize;
    maxSize   = _readMeanSize;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankSplitter::Iterator::Iterator(const BankSplitter& bank)
    : _dataRef (0), _itRef(0), _readMeanSize (bank._readMeanSize),
      _rank(0), _nbMax(0), _overlap(bank._overlap), _isDone(true)
{
    assert (bank._readMeanSize > 0);
    assert (bank._readMeanSize > bank._overlap);

    /** We get the first sequence of the referred bank. */
    setItRef (bank._reference->iterator());
    _itRef->first();
    assert (_itRef->isDone() == false);

    /** We create the reference Data object (that references the provided string).
     * NOTE : we force the encoding to be the same as the referred bank. */
    setDataRef (new Data (_itRef->item().getDataEncoding()));
    *_dataRef = (*_itRef)->getData();

    _offsetMax = _dataRef->size() - _readMeanSize;

    DEBUG (("refSize=%d  _readMeanSize=%d  _overlap=%d  _offsetMax=%d\n", _reference->size(), _readMeanSize, _overlap, _offsetMax));

    _offsets.clear();
    size_t idx;
    size_t delta = _readMeanSize -_overlap;
    for (idx=0; idx<_offsetMax; idx+=delta)
    {
        _offsets.push_back (make_pair (idx,_readMeanSize));
    }
    _offsets.push_back (make_pair (idx, _dataRef->size()-idx));

    DEBUG (("FOUND %d offsets\n", _offsets.size()));
    for (size_t i=0; i<_offsets.size(); i++)
    {
        DEBUG (("   %d  %d\n", _offsets[i].first, _offsets[i].second));
    }

    _nbMax = _offsets.size() * bank._coverage;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankSplitter::Iterator::~Iterator()
{
    setDataRef(0);
    setItRef(0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankSplitter::Iterator::first()
{
    _rank = -1;
    next ();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankSplitter::Iterator::next()
{
    _isDone = (++_rank >= _nbMax);
    if (!_isDone)
    {
        size_t offset = _offsets[_rank % (_offsets.size())].first;
        size_t size   = _offsets[_rank % (_offsets.size())].second;

        _item->getData().setRef (_dataRef, offset, size);
    }
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
