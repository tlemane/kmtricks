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

#include <gatb/kmer/impl/BloomAlgorithm.hpp>

#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/Banks.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/BloomBuilder.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/BagCache.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>

#include <iostream>
#include <map>
#include <math.h>

#include <gatb/tools/math/Integer.hpp>
#include <gatb/tools/math/NativeInt8.hpp>

#include <gatb/tools/storage/impl/StorageTools.hpp>

// We use the required packages
using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::math;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::storage::impl;

using namespace gatb::core::tools::math;

#define DEBUG(a)  //printf a

/********************************************************************************/
namespace gatb  {  namespace core  {   namespace kmer  {   namespace impl {
/********************************************************************************/

static const char* progressFormat1 = "Bloom: read solid kmers                ";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
BloomAlgorithm<span>::BloomAlgorithm (
    Storage& storage,
    Iterable<Count>*    solidIterable,
    size_t              kmerSize,
    float               nbitsPerKmer,
    size_t              nb_cores,
    BloomKind  bloomKind,
    IProperties*        options
)
    :  Algorithm("bloom", nb_cores, options),
       _kmerSize(kmerSize), _nbitsPerKmer(nbitsPerKmer), _bloomKind(bloomKind), _storage(storage), _solidIterable(0)
{
    setSolidIterable (solidIterable);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
BloomAlgorithm<span>::BloomAlgorithm (tools::storage::impl::Storage& storage)
    :  Algorithm("bloom", 0, 0),
       _kmerSize(0), _nbitsPerKmer(0), _storage(storage), _solidIterable(0)
{
    /** We get the kind in the storage. */
    string kind = _storage(this->getName()).getProperty ("kind");

    /** We parse the type. */
    parse (kind, _bloomKind);

    string xmlString = _storage(this->getName()).getProperty ("xml");
    stringstream ss; ss << xmlString;   getInfo()->readXML (ss);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
BloomAlgorithm<span>::~BloomAlgorithm ()
{
    setSolidIterable (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void BloomAlgorithm<span>::execute ()
{

    /** We get the number of solid kmers. */
    u_int64_t solidKmersNb = _solidIterable->getNbItems();

    float     NBITS_PER_KMER     = _nbitsPerKmer;
    u_int64_t estimatedBloomSize = (u_int64_t) (solidKmersNb * NBITS_PER_KMER);
    size_t    nbHash             = (int)floorf (0.7*NBITS_PER_KMER);

    if (estimatedBloomSize ==0 ) { estimatedBloomSize = 1000; }

    /** We create the kmers iterator from the solid file. */
    Iterator <Count>* itKmers = createIterator<Count> (
        _solidIterable->iterator(),
        solidKmersNb,
        progressFormat1
    );
    LOCAL (itKmers);

    /** We use a bloom builder. */
    BloomBuilder<span> builder (estimatedBloomSize, nbHash, _kmerSize, _bloomKind, getDispatcher()->getExecutionUnitsNumber());

    /** We instantiate the bloom object. */
    IBloom<Type>* bloom = 0;
    {
        TIME_INFO (getTimeInfo(), "build_from_kmers");
        bloom = builder.build (itKmers);
    }
    LOCAL (bloom);

    /** We save the bloom. */
    StorageTools::singleton().saveBloom<Type> (_storage.getGroup(this->getName()), "bloom", bloom, _kmerSize);

    /** We gather some statistics. */
    getInfo()->add (1, "stats");
    getInfo()->add (2, "kind",           "%s",  toString(_bloomKind));
    getInfo()->add (2, "bitsize",        "%ld", bloom->getBitSize());
    getInfo()->add (2, "nb_hash",        "%d",  bloom->getNbHash());
    getInfo()->add (2, "nbits_per_kmer", "%f",  _nbitsPerKmer);
    getInfo()->add (1, getTimeInfo().getProperties("time"));

    /** We save the kind in the storage. */
    _storage.getGroup(this->getName()).addProperty ("kind", toString(_bloomKind));
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
