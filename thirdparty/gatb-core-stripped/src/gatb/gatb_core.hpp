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

/** \file gatb_core.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Global header file
 */

#ifndef _GATB_CORE_HPP_
#define _GATB_CORE_HPP_

/********************************************************************************/

#include <gatb/system/impl/System.hpp>

#include <gatb/tools/designpattern/impl/Observer.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <gatb/tools/collections/impl/BagCache.hpp>
#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/BagPartition.hpp>
#include <gatb/tools/collections/impl/Bloom.hpp>
#include <gatb/tools/collections/impl/BloomGroup.hpp>
#include <gatb/tools/collections/impl/ContainerSet.hpp>
#include <gatb/tools/collections/impl/Hash16.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>
#include <gatb/tools/collections/impl/OAHash.hpp>
#include <gatb/tools/collections/impl/IterableHelpers.hpp>

#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/storage/impl/StorageTools.hpp>

#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/tools/misc/impl/Histogram.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/Tool.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>
#include <gatb/tools/misc/impl/StringLine.hpp>
#include <gatb/tools/misc/impl/LibraryInfo.hpp>
#include <gatb/tools/misc/impl/Tokenizer.hpp>

#include <gatb/tools/math/NativeInt64.hpp>
#include <gatb/tools/math/NativeInt128.hpp>
#include <gatb/tools/math/LargeInt.hpp>

#include <gatb/bank/impl/Banks.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>
#include <gatb/bank/impl/BankConverterAlgorithm.hpp>

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/BloomBuilder.hpp>
#include <gatb/kmer/impl/SortingCountAlgorithm.hpp>
#include <gatb/kmer/impl/BankKmers.hpp>
#include <gatb/kmer/impl/CountProcessor.hpp>

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::storage;
using namespace gatb::core::tools::storage::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::math;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

/********************************************************************************/

#endif /* _GATB_CORE_HPP_ */
