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

/** \file LinearCounter.hpp
 *  \brief Linear counter for kmers using a bloom
 */

#ifndef _GATB_CORE_KMER_IMPL_BLOOM_BUILDER_HPP_
#define _GATB_CORE_KMER_IMPL_BLOOM_BUILDER_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>

#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/collections/impl/Bloom.hpp>

#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/tools/misc/api/IProperty.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

#include <gatb/tools/math/NativeInt8.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

/*
 * This class is a linear counter
 */
template<size_t span=KMER_DEFAULT_SPAN> class LinearCounter
{
public:

    /** Shortcuts. */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::Count Count;
    typedef typename Kmer<span>::ModelCanonical::Kmer  Kmer;

    /** */
    LinearCounter (u_int64_t   bloom_size);

    void add (const Type& kmer);
    long count();
    bool is_accurate();

    ~LinearCounter();

private:

    gatb::core::tools::collections::impl::IBloom<Type>* bloom;

    u_int64_t _bloomSize;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMPL_LINEARCOUNTER_HPP_ */
