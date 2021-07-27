/*****************************************************************************
 *   kmtricks
 *   Authors: T. Lemane
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

#pragma once
#include <gatb/gatb_core.hpp>
#include <kmtricks/kmer.hpp>

namespace km {

template<size_t span>
inline void copy_gatb_kmers(Kmer<span>& kmtricks, const typename ::Kmer<span>::Type& gatb)
{
  for (size_t i=0; i<kmtricks.m_n_data; i++)
    kmtricks.get_data64_unsafe()[i] = gatb.get_data()[i];
};

template<>
inline void copy_gatb_kmers(Kmer<32>& kmtricks, const typename ::Kmer<32>::Type& gatb)
{
  kmtricks.set64(gatb.getVal());
}

#ifdef __SIZEOF_INT128__
template<>
inline void copy_gatb_kmers(Kmer<64>& kmtricks, const typename ::Kmer<64>::Type& gatb)
{
  kmtricks.set128(gatb.get_128());
}
#endif

inline void dump_pinfo(PartiInfo<5>* pinfo, uint32_t nb_parts, const std::string& path)
{
  std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
  for (uint32_t i=0; i<nb_parts; i++)
    out << std::to_string(pinfo->getNbKmer(i)) << "\n";
}

using props_t = std::shared_ptr<IProperties>;

inline props_t get_properties()
{
  props_t props = std::make_shared<Properties>();
  return props;
}

inline IProperties* get_config_properties(uint32_t kmer_size,
                                     uint32_t minim_size,
                                     uint32_t minim_type,
                                     uint32_t repart_type,
                                     uint32_t abundance_min,
                                     uint32_t nb_parts,
                                     uint32_t max_memory = 8000)
{
  //props_t props = std::make_shared<Properties>();
  IProperties* props = new Properties();
  props->add(0, "-kmer-size", "%d", kmer_size);
  props->add(0, "-minimizer-size", "%d", minim_size);
  props->add(0, "-minimizer-type", "%d", minim_type);
  props->add(0, "-repartition-type", "%d", repart_type);
  props->add(0, "-abundance-min", "%d", abundance_min);
  props->add(0, "-abundance-max", "%d", DMAX_C);
  props->add(0, "-solidity-kind", "sum");
  props->add(0, "-max-disk", "%d", 0);
  props->add(0, "-max-memory", "%d", max_memory);
  props->add(0, "-nb-cores", "%d", 1);
  props->add(0, "-storage-type", "0");
  props->add(0, "nb_partitions", "%d", nb_parts);
  return props;
}

inline props_t get_repart_properties()
{
  props_t props = std::make_shared<Properties>();
  return props;
}

inline props_t get_superk_properties()
{
  props_t props = std::make_shared<Properties>();
  return props;
}

inline props_t get_count_properties()
{
  props_t props = std::make_shared<Properties>();

  return props;
}

};