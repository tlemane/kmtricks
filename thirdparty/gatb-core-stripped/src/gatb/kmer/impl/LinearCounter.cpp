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

#include "LinearCounter.hpp"
#include <cmath> // for log2f

using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;


using namespace gatb::core::tools::math;

using namespace gatb::core::kmer::impl;


/********************************************************************************/
namespace gatb          {
namespace core          {
namespace kmer          {
namespace impl          {
/********************************************************************************/

    template<size_t span>
    LinearCounter<span>:: LinearCounter (
        u_int64_t   bloom_size
    )
        : _bloomSize (bloom_size) 
    {
        tools::misc::BloomKind bloomKind = tools::misc::BLOOM_BASIC;
        bloom = tools::collections::impl::BloomFactory::singleton().createBloom<Type> (
                                bloomKind, (u_int64_t)(bloom_size), 1, 31
                                                        );
    }

    template<size_t span>
    void LinearCounter<span>::add(const Type &kmer)
    {
        bloom->insert(kmer);
    }

    template<size_t span>
    long LinearCounter<span>::count()
    {
        long weight = bloom->weight();
        long count =  ( (-1.0*_bloomSize) * logf( (1.0*_bloomSize - weight) / _bloomSize ) );  // linear counter cardinality estimation

        bool debug=false;
        if (debug)
        {
            printf("linear counter load factor: %0.2f\n",(1.0*weight/_bloomSize));
            printf("weight: %ld bloomsize: %ld count: %ld\n",weight,(long)_bloomSize, count);
        }
        return count;
    }

    template<size_t span>
    bool LinearCounter<span>::is_accurate()
    {
        long weight = bloom->weight();
        float load_factor = (1.0*weight/_bloomSize);
        return load_factor < 0.99;
    }

    template<size_t span>
    LinearCounter<span>::~LinearCounter()
    {
        delete bloom;
    }

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/


