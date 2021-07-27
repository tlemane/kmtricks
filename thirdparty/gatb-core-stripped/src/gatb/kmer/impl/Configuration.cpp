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

#include <gatb/kmer/impl/Configuration.hpp>
#include <gatb/system/api/IMemory.hpp>

/********************************************************************************/

using namespace std;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace kmer          {
namespace impl          {
/********************************************************************************/

#define DEBUG(a)  //printf a

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Properties Configuration::getProperties() const
{
    Properties result ("config");

    result.add (1, "kmer_size",         "%ld", _kmerSize);
    result.add (1, "mini_size",         "%ld", _minim_size);
    result.add (1, "solidity_kind",      "%s", toString(_solidityKind).c_str());

	if(_solidityKind==KMER_SOLIDITY_CUSTOM)
	{
		std::stringstream ss2;
		for (size_t i=0; i<_solidVecUserNb; i++)  {  ss2 <<  " " << _solidVec[i]; }
		result.add (1, "custom_solidity",     ss2.str());
	}

	
    std::stringstream ss;
    for (size_t i=0; i<_abundanceUserNb; i++)  {  ss << (i==0 ? "" : " ") << _abundance[i].getBegin(); }

    result.add (1, "abundance_min",     ss.str());
    result.add (1, "abundance_max",     "%ld", _abundance[0].getEnd());

	

	
    result.add (1, "available_space",   "%ld", _available_space);
    result.add (1, "estimated_sequence_number",   "%ld", _estimateSeqNb);
    result.add (1, "estimated_sequence_volume",   "%ld", _estimateSeqTotalSize / system::MBYTE);
    result.add (1, "estimated_kmers_number",      "%ld", _kmersNb);
    result.add (1, "estimated_kmers_volume",      "%ld", _volume);
    result.add (1, "max_disk_space",    "%ld", _max_disk_space);
    result.add (1, "max_memory",        "%ld", _max_memory);
    result.add (1, "nb_passes",         "%d",  _nb_passes);
    result.add (1, "nb_partitions",     "%d",  _nb_partitions);
    result.add (1, "nb_bits_per_kmer",  "%d",  _nb_bits_per_kmer);
    result.add (1, "nb_cores",          "%d",  _nbCores);
    result.add (1, "minimizer_type",    "%s",  (_minimizerType == 0) ? "lexicographic (kmc2 heuristic)" : "frequency");
    result.add (1, "repartition_type",  "%s",  (_repartitionType == 0) ? "unordered" : "ordered");

    result.add (1, "nb_cores_per_partition",     "%d",  _nbCores_per_partition);
    result.add (1, "nb_partitions_in_parallel",  "%d",  _nb_partitions_in_parallel);
    result.add (1, "nb_cached_items_per_core_per_part",  "%d",  _nb_cached_items_per_core_per_part);

    result.add (1, "nb_banks",  "%d",  _nb_banks);

    return result;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Configuration::load (tools::storage::impl::Group& group)
{

	_isComputed = true;

    tools::storage::impl::Storage::istream is (group, "config");

    is.read ((char*)&_kmerSize,                    sizeof(_kmerSize));
    is.read ((char*)&_minim_size,                    sizeof(_minim_size));
    is.read ((char*)&_repartitionType,                    sizeof(_repartitionType));
    is.read ((char*)&_minimizerType,                    sizeof(_minimizerType));
    is.read ((char*)&_max_disk_space,                    sizeof(_max_disk_space));
    is.read ((char*)&_max_memory,                    sizeof(_max_memory));
    is.read ((char*)&_nbCores,                    sizeof(_nbCores));
    is.read ((char*)&_nb_partitions_in_parallel,                    sizeof(_nb_partitions_in_parallel));
    is.read ((char*)&_abundanceUserNb,                    sizeof(_abundanceUserNb));
    _abundance.resize (_abundanceUserNb);
    //is.read ((char*)_abundance.data(),    sizeof(tools::misc::CountRange)*_abundance.size());

    is.read ((char*)&_nbCores_per_partition,                sizeof(_nbCores_per_partition));
    is.read ((char*)&_estimateSeqNb,                    sizeof(_estimateSeqNb));
    is.read ((char*)&_estimateSeqTotalSize,             sizeof(_estimateSeqTotalSize));
    is.read ((char*)&_estimateSeqMaxSize,                sizeof(_estimateSeqMaxSize));
    is.read ((char*)&_available_space,    sizeof(_available_space));
    is.read ((char*)&_volume, sizeof(_volume));
    is.read ((char*)&_kmersNb,           sizeof(_kmersNb));
    is.read ((char*)&_nb_passes,           sizeof(_nb_passes));
    is.read ((char*)&_nb_partitions,           sizeof(_nb_partitions));
    is.read ((char*)&_nb_bits_per_kmer,           sizeof(_nb_bits_per_kmer));
    is.read ((char*)&_nb_banks,           sizeof(_nb_banks));
    is.read ((char*)&_nb_cached_items_per_core_per_part,           sizeof(_nb_cached_items_per_core_per_part));


}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Configuration::save (tools::storage::impl::Group& group)
{
    DEBUG (("[Config::save]\n"));

    tools::storage::impl::Storage::ostream os (group, "config");

    os.write ((const char*)&_kmerSize,                    sizeof(_kmerSize));
    os.write ((const char*)&_minim_size,                    sizeof(_minim_size));
    os.write ((const char*)&_repartitionType,                    sizeof(_repartitionType));
    os.write ((const char*)&_minimizerType,                    sizeof(_minimizerType));
    os.write ((const char*)&_max_disk_space,                    sizeof(_max_disk_space));
    os.write ((const char*)&_max_memory,                    sizeof(_max_memory));
    os.write ((const char*)&_nbCores,                    sizeof(_nbCores));
    os.write ((const char*)&_nb_partitions_in_parallel,                    sizeof(_nb_partitions_in_parallel));
    os.write ((const char*)&_abundanceUserNb,                    sizeof(_abundanceUserNb));
    //os.write ((const char*)_abundance.data(),    sizeof(tools::misc::CountRange)*_abundance.size());


    os.write ((const char*)&_nbCores_per_partition,                sizeof(_nbCores_per_partition));
    os.write ((const char*)&_estimateSeqNb,                    sizeof(_estimateSeqNb));
    os.write ((const char*)&_estimateSeqTotalSize,             sizeof(_estimateSeqTotalSize));
    os.write ((const char*)&_estimateSeqMaxSize,                sizeof(_estimateSeqMaxSize));
    os.write ((const char*)&_available_space,    sizeof(_available_space));
    os.write ((const char*)&_volume, sizeof(_volume));
    os.write ((const char*)&_kmersNb,           sizeof(_kmersNb));
    os.write ((const char*)&_nb_passes,           sizeof(_nb_passes));
    os.write ((const char*)&_nb_partitions,           sizeof(_nb_partitions));
    os.write ((const char*)&_nb_bits_per_kmer,           sizeof(_nb_bits_per_kmer));
    os.write ((const char*)&_nb_banks,           sizeof(_nb_banks));
    os.write ((const char*)&_nb_cached_items_per_core_per_part,           sizeof(_nb_cached_items_per_core_per_part));

    os.flush();

}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

