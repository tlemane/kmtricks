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

#ifndef _CONFIGURATION_HPP_
#define _CONFIGURATION_HPP_

/********************************************************************************/

#include <gatb/system/api/types.hpp>
#include <gatb/tools/misc/api/Enums.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/api/Range.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

class Configuration
{
public:

    /** */
    Configuration ()
    : _kmerSize(0), _minim_size(0), _repartitionType(0), _minimizerType(0),
      _solidityKind(tools::misc::KMER_SOLIDITY_SUM),
      _max_disk_space(0), _max_memory(0),
      _nbCores(0), _nb_partitions_in_parallel(0), _abundanceUserNb(0), _storage_type(tools::storage::impl::STORAGE_FILE) ,
      _isComputed(false), _nbCores_per_partition(0),
      _estimateSeqNb(0), _estimateSeqTotalSize(0), _estimateSeqMaxSize(0),
      _available_space(0), _volume(0), _kmersNb(0), _nb_passes(0), _nb_partitions(0), _nb_bits_per_kmer(0), _nb_banks(0) {}

    /****************************************/
    /**             PROVIDED                */
    /****************************************/

    size_t      _kmerSize;
    size_t      _minim_size;
    size_t      _repartitionType;
    size_t      _minimizerType;

    tools::misc::KmerSolidityKind _solidityKind;

    u_int64_t   _max_disk_space;
    u_int32_t   _max_memory;

    size_t      _nbCores;
    size_t      _nb_partitions_in_parallel;

    std::vector<tools::misc::CountRange>  _abundance;
    size_t _abundanceUserNb;

    tools::storage::impl::StorageMode_e _storage_type;
	
	std::vector<bool> _solidVec;
	size_t _solidVecUserNb;

    /****************************************/
    /**             COMPUTED                */
    /****************************************/

    bool        _isComputed;

    size_t      _nbCores_per_partition;

    u_int64_t   _estimateSeqNb;
    u_int64_t   _estimateSeqTotalSize;
    u_int64_t   _estimateSeqMaxSize;

    u_int64_t   _available_space;
    u_int64_t   _volume;
    u_int64_t   _kmersNb;

    u_int32_t   _nb_passes;
    u_int32_t   _nb_partitions;

    u_int16_t   _nb_bits_per_kmer;

    u_int16_t   _nb_banks;
    
    u_int32_t   _nb_cached_items_per_core_per_part;


    /****************************************/
    /**               MISC                  */
    /****************************************/
    tools::misc::impl::Properties getProperties() const;

    /** Load config properties from a storage object.
     * \param[in] group : group where the repartition table has to be loaded */
    void load (tools::storage::impl::Group& group);

    /** Save config properties into a storage object.
     * \param[in] group : group where the repartition table has to be saved */
    void save (tools::storage::impl::Group& group);
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _CONFIGURATION_HPP_ */
