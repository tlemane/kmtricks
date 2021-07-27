/*****************************************************************************
 *   kmtricks
 *   Authors: T. Lemane, R. Chikhi, GATB Team
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
#define NONCANONICAL
#include <gatb/gatb_core.hpp>

#include <gatb/kmer/impl/Sequence2SuperKmer.hpp>
#include <kmtricks/io/superk_storage.hpp>
#include <kmtricks/kmer.hpp>
namespace km {

template<size_t span>
class KmFillPartitions : public gatb::core::kmer::impl::Sequence2SuperKmer<span>
{
public:
  typedef typename Sequence2SuperKmer<span>::Type   Type;
  typedef typename Sequence2SuperKmer<span>::Model  Model;
  typedef typename Model::Kmer                      KmerType;
  typedef typename ::Kmer<span>::SuperKmer          SuperKmer;

  KmFillPartitions (Model& model,
                    size_t p,
                    size_t cp,
                    size_t nb_partitions,
                    size_t cache_items,
                    IteratorListener* progress,
                    BankStats& bank_stats,
                    Partition<Type>* partition,
                    Repartitor& repartition,
                    PartiInfo<5>& pinfo,
                    SuperKStorageWriter* superk)
    : Sequence2SuperKmer<span>(model, p, cp, nb_partitions, progress, bank_stats),
      m_kx(4),
      m_extern_pinfo(pinfo),
      m_local_pinfo(nb_partitions, model.getMmersModel().getKmerSize()),
      m_repartition(repartition),
      m_superk_files(superk)
  {
    m_mask_radix.setVal(static_cast<uint64_t>(255));
    m_mask_radix = m_mask_radix << ((this->_kmersize - 4) * 2);
  }

  void processSuperkmer(SuperKmer& superKmer)
  {
    if ((superKmer.minimizer % this->_nbPass) == this->_pass && superKmer.isValid())
    {
      size_t p = m_repartition(superKmer.minimizer);
      superKmer.save(p, m_superk_files);
      m_local_pinfo.incSuperKmer_per_minimBin(superKmer.minimizer, superKmer.size());

      Type radix_kxmer_forward, radix_kxmer;
      bool prev_which = superKmer[0].which();
      size_t kx_size = 0;
      radix_kxmer_forward = getHeavyWeight(superKmer[0].value());

      for (size_t ii=1; ii<superKmer.size(); ii++)
      {
        if (superKmer[ii].which() != prev_which || kx_size >= m_kx)
        {
          if (prev_which)
          {
            radix_kxmer = radix_kxmer_forward;
          }
          else
          {
            radix_kxmer = getHeavyWeight(superKmer[ii-1].value());
          }
          m_local_pinfo.incKmer_and_rad(p, radix_kxmer.getVal(), kx_size);
          radix_kxmer_forward = getHeavyWeight(superKmer[ii].value());
          kx_size = 0;
        }
        else
        {
          kx_size++;
        }
        prev_which = superKmer[ii].which();
      }
      if (prev_which)
      {
        radix_kxmer = radix_kxmer_forward;
      }
      else
      {
        radix_kxmer = getHeavyWeight(superKmer[superKmer.size()-1].value());
      }
      m_local_pinfo.incKmer_and_rad(p, radix_kxmer.getVal(), kx_size);
      this->_nbWrittenKmers += superKmer.size();
    }
  }

  virtual ~KmFillPartitions()
  {
    m_extern_pinfo.add_sync(m_local_pinfo);
  }

private:
  Type getHeavyWeight (const Type& kmer) const
  {
    return (kmer & m_mask_radix) >> ((this->_kmersize - 4)*2);
  }

private:
  size_t m_kx;
  PartiInfo<5>& m_extern_pinfo;
  PartiInfo<5>  m_local_pinfo;
  Type m_mask_radix;
  Repartitor& m_repartition;
  SuperKStorageWriter* m_superk_files;
};


}