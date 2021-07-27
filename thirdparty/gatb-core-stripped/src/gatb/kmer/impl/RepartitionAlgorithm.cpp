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

#include <gatb/kmer/impl/RepartitionAlgorithm.hpp>

#include <gatb/system/impl/System.hpp>

#include <gatb/bank/impl/Banks.hpp>
#include <gatb/bank/impl/BankHelpers.hpp>

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/BloomBuilder.hpp>
#include <gatb/kmer/impl/BankKmers.hpp>
#include <gatb/kmer/impl/Sequence2SuperKmer.hpp>

#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/BagCache.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>

#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>

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

static const char* progressFormat0 = "DSK: Collecting stats on %s ";


/********************************************************************************/

template<size_t span>
class MmersFrequency
{
public:
    /** Shortcut. */
    typedef typename RepartitorAlgorithm<span>::ModelCanonical ModelCanonical;
    typedef typename ModelCanonical::Kmer                       KmerTypeCanonical;
    typedef typename RepartitorAlgorithm<span>::ModelDirect    ModelDirect;
    typedef typename ModelDirect::Kmer                       KmerTypeDirect;

    void operator() (Sequence& sequence)
    {
        /** We first check whether we got mmers from the sequence or not. */
        if (_minimodel.build (sequence.getData(), _mmers) == false)  { return; }

        /** We loop over the mmers of the sequence. */
        for (size_t i=0; i<_mmers.size(); i++)
        {
            if (_mmers[i].isValid() == false)
                continue;

            /** increment m-mer count */
            _m_mer_counts[_mmers[i].value().getVal()] ++;
        }

        if (_nbProcessedMmers > 500000)   {  _progress.inc (_nbProcessedMmers);  _nbProcessedMmers = 0;  }

        // see if we need to stop. 
        // mimics SampleRepart below. but actually, as TruncatedIterator would have worked too, since all seqs are at least larger than a mmer. oh well. using CancellableIterator anyway.
        _nbSeqsSeenSoFar ++;
        if (_nbSeqsSeenSoFar > _nbSeqsToSee)
        {
            *_cancelIterator = true;
        }

    }

    /** Constructor. */
    MmersFrequency (int mmerSize,  IteratorListener* progress,  uint32_t* m_mer_counts, size_t nbSeqsToSee /* for when to stop estimation*/, bool* cancelIterator)
    :
        _minimodel(mmerSize), _progress (progress,System::thread().newSynchronizer()),
      _m_mer_counts(m_mer_counts),  _nbProcessedMmers(0), _nbSeqsToSee(nbSeqsToSee), _nbSeqsSeenSoFar(0), _cancelIterator(cancelIterator)
    {
        u_int64_t nbminim = (uint64_t)pow(4.0,mmerSize);

        for (u_int64_t i = 0; i < nbminim; i++)  {   _m_mer_counts[i] = 0;  }
    }

protected:

    ModelCanonical           _minimodel;
    //ModelDirect                _minimodel;
    vector<KmerTypeCanonical>  _mmers;
    //vector<KmerTypeDirect>  _mmers;
    ProgressSynchro         _progress;
    uint32_t*               _m_mer_counts;
    size_t                  _nbProcessedMmers;
    size_t        _nbSeqsToSee;
    size_t        _nbSeqsSeenSoFar;
    bool *            _cancelIterator;
};

/********************************************************************************/
/* This functor class takes a Sequence as input, splits it into super kmers and
 * get information about the distribution of minimizers.
 */
template<size_t span>
class SampleRepart  : public Sequence2SuperKmer<span>
{
    //ie ce sera posible d avoir plus dinfo , estim ram max par exemple ?
public:

    /** Shortcut. */
    typedef typename Sequence2SuperKmer<span>::Type             Type;
    typedef typename Sequence2SuperKmer<span>::Model            Model;
    typedef typename Model::Kmer                           KmerType;
    typedef typename Kmer<span>::SuperKmer                 SuperKmer;

    /** */
    void processSuperkmer (SuperKmer& superKmer)
    {
        DEBUG (("SampleRepart: should count superk %i \n", superKmer.size()));

        if ((superKmer.minimizer % this->_nbPass) == this->_pass && superKmer.isValid() ) //check if falls into pass
        {
            bool prev_which = superKmer[0].which();
            size_t kx_size = 0;

            /** Shortcut. */
            size_t superKmerLen = superKmer.size();

            /** We increase superkmer counter the current minimizer. */
            _pInfo.incSuperKmer_per_minimBin (superKmer.minimizer, superKmerLen);

            /** We loop over the kmer of the superkmer (except the first one).
             *  We update the pInfo each time we find a kxmer in the superkmer. */
            for (size_t ii=1 ; ii < superKmerLen; ii++)
            {
                /** A kxmer is defined by having successive canonical kmers. Here, we just care that
                 * successive kmer values are on the same strand. */
                if (superKmer[ii].which() != prev_which || kx_size >= _kx) // kxmer_size = 1 //cost should diminish with larger kxmer
                {
                    /** We increase the number of kxmer found for the current minimizer. */
                    _pInfo.incKxmer_per_minimBin (superKmer.minimizer);
                    kx_size = 0;
                }
                else
                {
                    kx_size++;
                }

                prev_which = superKmer[ii].which() ;
            }

            /** We add the pending kxmer to the bin. */
            _pInfo.incKxmer_per_minimBin (superKmer.minimizer);
        
            // see if we need to stop
            _nbSuperKmersSeenSoFar ++;
            if (_nbSuperKmersSeenSoFar > _nbSeqsToSee) // it's a bit of an approximation here (comparing # of superkmers to # of seqs)
            {
                //cout << _nbSuperKmersSeenSoFar << "  " <<  _nbSeqsToSee << endl;
                *_cancelIterator = true;
            }
        }
    }

    /** Constructor. */
    SampleRepart (
        Model&            model,
        Configuration&    config,
        size_t            nbPartitions,
        IteratorListener* progress,
        bool *            cancelIterator,
        size_t            nbSeqsToSee,
        BankStats&        bankStats,
        PartiInfo<5>&     pInfo
    )
    :   Sequence2SuperKmer<span> (model, 1, 0, nbPartitions, progress, bankStats)
        ,_kx(4), _pInfo(pInfo),
        _cancelIterator(cancelIterator), _nbSeqsToSee(nbSeqsToSee), _nbSuperKmersSeenSoFar(0)
    {
    }


private:
    size_t        _kx;
    PartiInfo<5>& _pInfo;

    bool*         _cancelIterator;
    size_t        _nbSeqsToSee;
    size_t        _nbSuperKmersSeenSoFar;
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
RepartitorAlgorithm<span>::RepartitorAlgorithm (
    IBank* bank,
    Group& group,
    const Configuration& config,
    unsigned int nb_cores,
    tools::misc::IProperties*   options
)
    :  Algorithm("repartition", nb_cores, options), _config(config), _bank(bank), _group(group), _freq_order(0)
{
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
RepartitorAlgorithm<span>::~RepartitorAlgorithm ()
{
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
void RepartitorAlgorithm<span>::execute ()
{
    /** We compute the distribution of the minimizers. As a result, we will have a hash function
     * that gives a hash code for a minimizer value.
     * IMPORTANT ! we have to give the passes number because it has impact on the computation. */
    Repartitor repartitor (_config._nb_partitions, _config._minim_size, _config._nb_passes);

    /* now is a good time to switch to frequency-based minimizers if required:
      because right after we'll start using minimizers to compute the distribution
      of superkmers in bins */
    if (_config._minimizerType == 1)  {  computeFrequencies (repartitor);  }

    computeRepartition (repartitor);
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
void RepartitorAlgorithm<span>::computeFrequencies (Repartitor& repartitor)
{
    DEBUG (("RepartitorAlgorithm<span>::computeFrequencies\n"));

    u_int64_t estimateSeqNb;
    u_int64_t estimateSeqTotalSize;
    u_int64_t estimateSeqMaxSize;

    _bank->estimate (estimateSeqNb, estimateSeqTotalSize, estimateSeqMaxSize);

    u_int64_t nbseq_sample = std::min ( u_int64_t (estimateSeqNb * 0.05) ,u_int64_t( 50000000ULL) ) ;
    // TODO would be better to just stop estimating minimizer frequency when it becomes stable. not after a fixed number of reads

    if (nbseq_sample == 0)
        nbseq_sample = 1;

    u_int64_t rg = ((u_int64_t)1 << (2*_config._minim_size));
    //cout << "\nAllocating " << ((rg*sizeof(uint32_t))/1024) << " KB for " << _minim_size <<"-mers frequency counting (" << rg << " elements total)" << endl;
    uint32_t *m_mer_counts = new uint32_t[rg];

    Model model (_config._kmerSize, _config._minim_size);

    Iterator<Sequence>* bank_it = _bank->iterator();
    LOCAL(bank_it);
    CancellableIterator<Sequence>* cancellable_it = new CancellableIterator<Sequence> (*bank_it);
    LOCAL(cancellable_it);

    /** We create a sequence iterator and give it a progress message */
    Iterator<Sequence>* it_all_reads = createIterator<Sequence> (
            cancellable_it,
            _bank->getNbItems(),
             "Approximating frequencies of minimizers"
            );
    LOCAL (it_all_reads);

    /** We compute an estimation of minimizers frequencies from a part of the bank. */
    
    SerialDispatcher serialDispatcher;
    serialDispatcher.iterate (it_all_reads,  MmersFrequency<span> (
        _config._minim_size, 0 /*_progress*/, m_mer_counts, 
        nbseq_sample,
        &(cancellable_it->_cancel))// will be set to true when iteration needs to be stopped
    );

    // single threaded, for debugging
    /*MmersFrequency<span> mmersfrequency(model, _progress, bstatsDummy, m_mer_counts);
    it_sample->iterate(mmersfrequency);*/

    /* sort frequencies */
    for (u_int64_t i(0); i < rg; i++)
    {
        if (m_mer_counts[i] > 0)
            _counts.push_back(make_pair(m_mer_counts[i],i));
    }
    delete[] m_mer_counts;

    sort (_counts.begin(), _counts.end());

    /* assign frequency to minimizers */
    _freq_order = new uint32_t[rg];

    for (u_int64_t i = 0; i < rg ; i++)
        _freq_order[i] = rg; // set everything not seen to highest value (not a minimizer)

    for (unsigned int i = 0; i < _counts.size(); i++)
    {
        _freq_order[_counts[i].second] = i;
    }

    // small but necessary trick: the largest minimizer has to have largest rank, as it's used as the default "largest" value
    _freq_order[rg-1] = rg-1;

    repartitor.setMinimizerFrequencies (_freq_order);
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
void RepartitorAlgorithm<span>::computeRepartition (Repartitor& repartitor)
{
    DEBUG (("RepartitorAlgorithm<span>::computeRepartition\n"));

    /** We create a kmer model; using the frequency order if we're in that mode */
    Model model (_config._kmerSize, _config._minim_size, typename Kmer<span>::ComparatorMinimizerFrequencyOrLex(), _freq_order);

    int mmsize = model.getMmersModel().getKmerSize();

    PartiInfo<5> sample_info (_config._nb_partitions, mmsize);

    string bankShortName = System::file().getBaseName(_bank->getId());

    // In case of multi bank counting, we get a sample from each bank
    if(_bank->getCompositionNb() > 1){

    	//cout << _bank->getCompositionNb() << endl;
    	SerialDispatcher serialDispatcher;
        BankStats bstatsDummy;

        u_int64_t nbseq_sample = (_config._estimateSeqNb / _config._nb_banks) * 0.01;
        nbseq_sample = max((u_int64_t)nbseq_sample, (u_int64_t)100000);

        Iterator<Sequence>* it = _bank->iterator(); LOCAL (it);
        std::vector<Iterator<Sequence>*> itBanks =  it->getComposition(); 

        for(size_t i=0; i<_config._nb_banks; i++){
            CancellableIterator<Sequence>* cancellable_it = new CancellableIterator<Sequence> (*itBanks[i]);
            LOCAL(cancellable_it);

    		/** We compute a distribution of Superkmers from a part of the bank. */
    		serialDispatcher.iterate (cancellable_it, SampleRepart<span> (
    			model,
    			_config,
    			_config._nb_partitions,
    			NULL,
    			&(cancellable_it->_cancel), // will be set to true when iteration needs to be stopped
    			nbseq_sample, // how many sequences we need to see
    			bstatsDummy,
    			sample_info
    		));

    		//cout << "end" << endl;
    		//it_all_reads->finalize() ;
    		//cancellable_it->finalize() ;
    		itBanks[i]->finalize();
    		//delete cancellable_it;
        }
    }
    else{

        Iterator<Sequence>* it = _bank->iterator();      LOCAL (it);
        CancellableIterator<Sequence>* cancellable_it = new CancellableIterator<Sequence> (*(it));
        LOCAL(cancellable_it);

		// how many seqs we need to see
		u_int64_t nbseq_sample = std::max ( u_int64_t (_config._estimateSeqNb * 0.05) ,u_int64_t( 1000000ULL) ) ;

		/** We create a sequence iterator and give it a progress message */
		Iterator<Sequence>* it_all_reads = createIterator<Sequence> (
				cancellable_it,
				_bank->getNbItems(),
				Stringify::format (progressFormat0, bankShortName.c_str()).c_str()
				);
		LOCAL (it_all_reads);

		BankStats bstatsDummy;

		/** We compute a distribution of Superkmers from a part of the bank. */
		SerialDispatcher serialDispatcher;
		serialDispatcher.iterate (it_all_reads, SampleRepart<span> (
			model,
			_config,
			_config._nb_partitions,
			NULL,
			&(cancellable_it->_cancel), // will be set to true when iteration needs to be stopped
			nbseq_sample, // how many sequences we need to see
			bstatsDummy,
			sample_info
		));
    }

    if (_config._minimizerType == 1)
    {
        repartitor.justGroup (sample_info, _counts);
    }
    else
    {
        repartitor.computeDistrib (sample_info);
        if (_config._repartitionType == 1)
        {
            repartitor.justGroupLexi (sample_info); // For bcalm, i need the minimizers to remain in order. so using this suboptimal but okay repartition
        }
    }

    /** We save the distribution (may be useful for debloom for instance). */
    repartitor.save (getGroup());
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
