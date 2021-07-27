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

/** \file BloomAlgorithm.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Bloom algorithm, ie. compute a Bloom filter from a set of reads
 */

#ifndef _BLOOM_ALGORITHM_HPP_
#define _BLOOM_ALGORITHM_HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Algorithm.hpp>

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/collections/impl/Bloom.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>

#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

template<size_t span=KMER_DEFAULT_SPAN>
class BloomAlgorithm : public gatb::core::tools::misc::impl::Algorithm
{
public:

    /** Shortcuts. */
    typedef typename kmer::impl::Kmer<span>::ModelCanonical Model;
    typedef typename kmer::impl::Kmer<span>::Type           Type;
    typedef typename kmer::impl::Kmer<span>::Count          Count;

    /** */
    BloomAlgorithm (
        tools::storage::impl::Storage&       storage,
        tools::collections::Iterable<Count>* solidIterable,
        size_t                               kmerSize,
        float                                nbitsPerKmer,
        size_t                               nb_cores = 0,
        tools::misc::BloomKind               bloomKind = tools::misc::BLOOM_DEFAULT,
        tools::misc::IProperties*            options    = 0
    );

    /** */
    BloomAlgorithm (tools::storage::impl::Storage& storage);

    /** */
    ~BloomAlgorithm ();

    /** */
    void execute ();

private:

    /** */
    size_t _kmerSize;

    /** */
    float _nbitsPerKmer;

    /** */
    tools::misc::BloomKind _bloomKind;

    /** */
    tools::storage::impl::Storage& _storage;

    /** */
    tools::collections::Iterable<Count>* _solidIterable;
    void setSolidIterable (tools::collections::Iterable<Count>* solidIterable)  {  SP_SETATTR(solidIterable); }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _BLOOM_ALGORITHM_HPP_ */

