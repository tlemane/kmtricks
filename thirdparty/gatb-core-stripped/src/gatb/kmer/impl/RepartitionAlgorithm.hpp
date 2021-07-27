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

/** \file RepartitionAlgorithm.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Repartition algorithm, ie. compute statistics on kmers
 */

#ifndef _REPARTITOR_ALGORITHM_HPP_
#define _REPARTITOR_ALGORITHM_HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/bank/api/IBank.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/Configuration.hpp>
#include <gatb/kmer/impl/PartiInfo.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>
#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

template<size_t span=KMER_DEFAULT_SPAN>
class RepartitorAlgorithm : public gatb::core::tools::misc::impl::Algorithm
{
public:

    typedef typename kmer::impl::Kmer<span>::ModelDirect                                ModelDirect;
    typedef typename kmer::impl::Kmer<span>::ModelCanonical                             ModelCanonical;
#ifdef NONCANONICAL
    typedef typename kmer::impl::Kmer<span>::template ModelMinimizer <ModelDirect>   Model;
#else
    typedef typename kmer::impl::Kmer<span>::template ModelMinimizer <ModelCanonical>   Model;
#endif

    /** */
    RepartitorAlgorithm (
        gatb::core::bank::IBank*        bank,
        tools::storage::impl::Group&    group,
        const Configuration&            config,
        unsigned int                    nb_cores = 0,
        tools::misc::IProperties*   options    = 0
    );

    /** */
    ~RepartitorAlgorithm ();

    /** */
    void execute ();

private:

    void computeFrequencies (Repartitor& repartitor);
    void computeRepartition (Repartitor& repartitor);

    Configuration _config;

    gatb::core::bank::IBank*      _bank;
    tools::storage::impl::Group&  _group;

    uint32_t* _freq_order;

    tools::storage::impl::Group& getGroup() { return  _group; }

    std::vector<std::pair<int, int> > _counts;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _REPARTITOR_ALGORITHM_HPP_ */

