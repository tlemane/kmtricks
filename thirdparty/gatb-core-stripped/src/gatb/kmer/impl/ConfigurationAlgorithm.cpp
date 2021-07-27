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

#include <gatb/kmer/impl/ConfigurationAlgorithm.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/collections/impl/OAHash.hpp>
#include <gatb/tools/misc/api/StringsRepository.hpp>
#include <gatb/tools/misc/impl/Tokenizer.hpp>
#include <gatb/kmer/impl/LinearCounter.hpp>

#include <cmath>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace kmer          {
namespace impl          {
/********************************************************************************/

using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

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

// estimated the number of distinct kmers in a dataset
// wrapper around a Linear Counter. Adapted from Kmergenie code.
// why not a hyperloglog? it seems that the transition from the 32bit-hash initial implementation to 64bits, and supporting billions of elements, is nontrivial, so i didn't bother
// probably deserves to be in its own file
template<size_t span>
class EstimateNbDistinctKmers
{
public:

    /** Shortcut. */
    typedef typename Kmer<span>::Type  Type;
    typedef typename Kmer<span>::ModelDirect     ModelDirect;
    typedef typename Kmer<span>::ModelCanonical  ModelCanonical;
#ifdef NONCANONICAL 
    typedef typename Kmer<span>::template ModelMinimizer <ModelDirect>   Model;
#else
    typedef typename Kmer<span>::template ModelMinimizer <ModelCanonical>   Model;
#endif
    typedef typename Model::Kmer                        KmerType;

    /** */
    void estimate()
    {
        nb_distinct_kmers =(unsigned long)( (float)(linearCounter->count( )) * ((float)nbKmersTotal / (float)nbProcessedKmers)); // dubious linear extrapolation, that's all I got

        abs_error = abs((long)(nb_distinct_kmers-previous_nb_distinct_kmers));

        previous_nb_distinct_kmers = nb_distinct_kmers;
    }

    /** */
    void operator() (Sequence& sequence)
    {
        /** We build the kmers from the current sequence. */
        if (model.build (sequence.getData(), kmers) == false)  {  throw "reached EOF"; return; }

        /** We loop over the kmers. */
        for (size_t i=0; i<kmers.size(); i++)
        {
            linearCounter->add((kmers[i].value()));


            // heuristics to stop early, i found that it's inaccurate with low coverage (e.g. on dsk/test/FiftyK.fastq)
            /*
            if (nbProcessedReads % eval_every_N_reads == 0 )
            {

                // let's see if the estimation converges..
                // the following stopping condition will grossly over-estimate the number of distinct kmers
                // but I expect the correct result to be in the same order of magnitude
                // and better to overestimate than underestimate (for both dsk and kmergenie)
                   estimate();
                   bool debug = true;
                   if (debug)
                       printf("linear estimator at %ld kmers, number of distinct kmers estimated now: %ld, abs error: %ld\n",nbProcessedKmers, nb_distinct_kmers, abs_error);
                   if (abs_error < previous_nb_distinct_kmers/20) // 5% error
                   {
                       throw "LinearCounter converged"; // well, "converged" is a big word
                       return;
                   }
                   if (!linearCounter->is_accurate())
                   {
                   printf("LinearCounter is inaccurate";
                   return;
                   }

            }*/

        }
        nbProcessedKmers += kmers.size();
        nbProcessedReads++;
        //if (nbProcessedReads % 100000 == 0) printf("nb: %ld\n",nbProcessedReads);

        // disabled progress
        //if (nbCurProgressKmers > 500000)   {  _progress.inc (nbCurProgressKmers);  nbCurProgressKmers = 0;  }
    }

    EstimateNbDistinctKmers (Model& model, u_int32_t max_memory, unsigned long nb_kmers_total, tools::dp::IteratorListener* progress)
        : model(model),  eval_every_N_reads(10000000),   nbKmersTotal(nb_kmers_total),
        nbProcessedKmers(0), nbCurProgressKmers(0), previous_nb_distinct_kmers(0), nbProcessedReads(0), abs_error(0)
        //, _progress  (progress,System::thread().newSynchronizer())
    {
        unsigned long size_linearCounter; // (in bits)
        /* let's set it to just use half of all memory available at most, ok? this isn't very robust for huge dataset, so to be tested*/
        /* if it's a tiny dataset, let's set it to total number of kmers */
        size_linearCounter = std::min( nb_kmers_total, (unsigned long) (max_memory*8*1024*1024/2) );
        linearCounter =  new LinearCounter<span>(size_linearCounter);
    }

    unsigned long getEstimation()
    {
        estimate();
        // soo.. if it's not accurate, let's assume we have a hugeload of kmers, and let's be safe, we return the total number of kmers
        if (!linearCounter->is_accurate())
        {
            cout << "Warning: linear counter was not accurate, returning worst-case estimation of number of distinct kmers";
            return nbKmersTotal;
        }
        return nb_distinct_kmers;
    }

private:

    /** Local resources. */
    Model&    model;
    unsigned long nbProcessedReads, nbProcessedKmers;
    unsigned long nbCurProgressKmers;
    unsigned long nbKmersTotal;
    unsigned long abs_error;
    vector<KmerType> kmers;
    LinearCounter<span> *linearCounter;
    int eval_every_N_reads;
    unsigned long previous_nb_distinct_kmers, nb_distinct_kmers;
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
ConfigurationAlgorithm<span>::ConfigurationAlgorithm (bank::IBank* bank, IProperties* input)
    : Algorithm("configuration", -1, input), _bank(0), _input (0)
{
    setBank  (bank);
    setInput (input);

    _config._kmerSize           = input->getInt (STR_KMER_SIZE);
    _config._minim_size         = input->getInt (STR_MINIMIZER_SIZE);
    _config._repartitionType    = input->getInt (STR_REPARTITION_TYPE);
    _config._minimizerType      = input->getInt (STR_MINIMIZER_TYPE);

    parse (input->getStr (STR_SOLIDITY_KIND), _config._solidityKind);

    _config._max_disk_space     = input->getInt (STR_MAX_DISK);
    _config._max_memory         = input->getInt (STR_MAX_MEMORY);
    _config._nbCores            = input->get(STR_NB_CORES) ? input->getInt(STR_NB_CORES) : 0;

    _config._abundance = getSolidityThresholds(input);
	
	if( _config._solidityKind == KMER_SOLIDITY_CUSTOM )
		_config._solidVec = getSolidityCustomVector(input);

    if (_config._nbCores == 0)  { _config._nbCores = system::impl::System::info().getNbCores(); }

    _config._nb_partitions_in_parallel = _config._nbCores;

    _config._nb_bits_per_kmer = Type::getSize();
    
    std::string storage_type = input->getStr(STR_STORAGE_TYPE);
    _config._storage_type =  tools::storage::impl::STORAGE_FILE;
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
ConfigurationAlgorithm<span>::~ConfigurationAlgorithm ()
{
    setBank  (0);
    setInput (0);
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
void ConfigurationAlgorithm<span>::execute ()
{
    /** By default, we want to have mmers of size 8. However (for unit tests for instance),
     * we may need to have kmer sizes less than 8; in such a case, we set by convention m=k-1. */
    if (_config._minim_size == 0)  {   _config._minim_size = 8;  }

    _config._minim_size = std::min ((int)_config._kmerSize-1, (int)_config._minim_size);

    /** We get some information about the bank. */
    _bank->estimate (_config._estimateSeqNb, _config._estimateSeqTotalSize, _config._estimateSeqMaxSize);
    
    //printf("_estimateSeqNb %llu _estimateSeqTotalSize %llu _estimateSeqMaxSize %llu \n", _config._estimateSeqNb, _config._estimateSeqTotalSize, _config._estimateSeqMaxSize);

    /** We get the number of sub banks. */
    _config._nb_banks = _bank->getCompositionNb();

	if(_config._nb_banks == 1 )
	{
		_config._solidityKind = KMER_SOLIDITY_SUM;
	}
	
    /** We memorize the number of abundance min values set by the user.
     * Note that it can be lower than the number of banks. */
    _config._abundanceUserNb = _config._abundance.size();

    if (_config._abundanceUserNb == 0)  {  throw system::Exception("Kmer solidity has no defined value");  }

    if (_config._abundanceUserNb > _config._nb_banks)
    {
        throw system::Exception ("Kmer solidity has more thresholds (%d) than banks (%d)",  _config._abundanceUserNb, _config._nb_banks);
    }

	_config._solidVecUserNb = _config._solidVec.size();

	if (_config._solidityKind == KMER_SOLIDITY_CUSTOM &&  _config._solidVecUserNb != _config._nb_banks)
	{
		throw system::Exception ("Kmer solidity custom has different number of values  (%d) than banks (%d)",  _config._solidVecUserNb, _config._nb_banks);
	}
	
	if (_config._solidityKind != KMER_SOLIDITY_CUSTOM )
	{
		_config._solidVec = std::vector<bool> (_config._nb_banks,true);
	}
	
    /** We complete missing thresholds with the value of the last one. */
    if (_config._abundanceUserNb < _config._nb_banks)
    {
        CountNumber lastValueMin = _config._abundance[_config._abundanceUserNb-1].getBegin();
        CountNumber lastValueMax = _config._abundance[_config._abundanceUserNb-1].getEnd();

        for (size_t i=_config._abundanceUserNb; i<_config._nb_banks; i++)
        {
            _config._abundance.push_back (tools::misc::CountRange (lastValueMin, lastValueMax));
        }
    }

    /** Some checks. */
    if (_config._estimateSeqNb==0)  { throw Exception ("Empty bank"); }

    // We get the available space (in MBytes) of the current directory.
    u_int64_t available_space_min = 2000;
    _config._available_space = System::file().getAvailableSpace (System::file().getCurrentDirectory()) / 1024;

    size_t meanSeqLen = (size_t) ( (double) _config._estimateSeqTotalSize / (double) _config._estimateSeqNb);
    size_t usedSeqLen = meanSeqLen > _config._kmerSize ? meanSeqLen : _config._kmerSize; // the latter used to be estimated max seq size, but that was too big of an overestimation in some cases (think one large sequence and plenty of sequences of length k-1). Let's see if setting to kmerSize works.

    int64_t kmersNb = (usedSeqLen - _config._kmerSize + 1) * _config._estimateSeqNb;

    /** We have to be sure that the kmers number is ok. */
    if (kmersNb <= 0)  {  throw Exception ("Configuration failed : estimated that longest sequence is %ld nt but kmer size is %ld", _config._estimateSeqMaxSize, _config._kmerSize);     }

    /** The estimated kmers number is ok. */
    _config._kmersNb  = kmersNb;

    _config._volume =  _config._kmersNb * sizeof(Type) / MBYTE;  // in MBytes

    //printf("estimated usedSeqLen %llu _estimateSeqNb %llu _kmersNb %llu  _volume %llu  \n", usedSeqLen, _config._estimateSeqNb, kmersNb, _config._volume);

    if (_config._volume == 0)   { _config._volume = 1; }    // tiny files fix

    u_int64_t volume_minim = _config._volume * 0.5 *1.2  ; //0.5 for using kxmers   1.2 if bad repartition of minimizers ( todo sampling to assert ram usage)

    if (volume_minim == 0)   { volume_minim = 1; }    // tiny files fix
    // volume_minim is used a bit later

    /** We get max(75%, 100% - X GB) */
    if (_config._max_disk_space == 0)  { _config._max_disk_space = std::max ((75*_config._available_space)/100, _config._available_space-available_space_min);  }
    if (_config._max_disk_space == 0)  { _config._max_disk_space = 10000; }

    if (_config._max_memory == 0)  {  _config._max_memory = System::info().getMemoryProject(); }
    if (_config._max_memory == 0)  {  _config._max_memory = 5000; }

    /* make sure to not use more mem than system, when max_memory has default value (useful for docker images) */
    if (_config._max_memory == 5000)  {
        unsigned long system_mem = System::info().getMemoryPhysicalTotal() / MBYTE;
        if (_config._max_memory > (system_mem * 2) / 3)
        {
            _config._max_memory = (system_mem * 2) / 3;
            if (_config._max_memory < 4500) // don't print a message if we were close to 5 GB anyway
                cout << "Warning: default memory usage (5000 MB) is close or above system max, setting memory to: " << _config._max_memory << " MB" << endl;
        }
    }

    assert (_config._max_disk_space > 0);

    _config._nb_passes = ( (_config._volume/4) / _config._max_disk_space ) + 1; //minim, approx volume /switched to approx /4 (was/3) because of more efficient superk storage
    //_nb_passes = 1; //do not constrain nb passes on disk space anymore (anyway with minim, not very big)
    //increase it only if ram issue

    //printf("_volume  %lli volume_minim %lli _max_disk_space %lli  _nb_passes init %i  \n", _config._volume,volume_minim, _config._max_disk_space, _config._nb_passes);
    size_t max_open_files = System::file().getMaxFilesNumber() / 2;


    if (_config._storage_type == tools::storage::impl::STORAGE_FILE)
    {
        //std::cout << "using less max_open_open files (" << max_open_files << "), by 3x, due to storage file setting" << std::endl;
        max_open_files /= 3; // will need to open twice in STORAGE_FILE instead of h d f 5, so this adjustment is needed. needs to be fixed later by putting partitions inside the same file. but i'd rather not do it in the current messy collection/group/partition h d f 5-inspired system. overall, that's a FIXME
    }

#if 0
    /* disabled by default; this was an experiment */
    float est_volume_distinct_ratio;
    if (_flagEstimateNbDistinctKmers)
    {
        /* we estimate the volume of distinct kmers vs total number of kmers.
         * we store it in the variable "est_volume_distinct_ratio"
         * to compute it, we need a linear counter, let's call it now */

        TIME_INFO (getTimeInfo(), "estimate_distinct_kmers");
        Iterator<Sequence>* itSeq = _bank->iterator();
        LOCAL (itSeq);

        //_progress->setMessage (progressFormat0); // not touching progress here anymore
        Model model (_kmerSize, _minim_size);
        EstimateNbDistinctKmers<span> estimate_nb_distinct_kmers_function(model, _max_memory, kmersNb, _progress);

        /** We launch the iteration of the sequences iterator with the created functors. */
        try {
            itSeq->iterate (estimate_nb_distinct_kmers_function);
        }
        catch (const char* except)
        {

        }
        _estimatedDistinctKmerNb = estimate_nb_distinct_kmers_function.getEstimation();
        est_volume_distinct_ratio = (float) _estimatedDistinctKmerNb / (float)kmersNb;
        //est_volume_distinct_ratio = 1; // for debug
        /* est_volume_distinct_ratio == 1 mean that we guarantee worst case the memory usage,
           the value mean that, on average, a k-mer will be seen 'est_volume_distinct_ratio' times */
        // if wrongly estimated, the error 'OAHash: max rehashes..' can happen
        printf ("LinearCounter done, estimated %ld number of distinct kmers, ratio to total number of kmers: %.2f\n", (long)_estimatedDistinctKmerNb, est_volume_distinct_ratio);
    }
#endif
    u_int64_t volume_per_pass;
    do  {

        assert (_config._nb_passes > 0);
        volume_per_pass = volume_minim / _config._nb_passes;

        assert (_config._max_memory > 0);
        //printf("volume_per_pass %lli  _nbCores %zu _max_memory %i \n",volume_per_pass, _nbCores,_max_memory);

        // _nb_partitions  = ( (volume_per_pass*_nbCores) / _max_memory ) + 1;
        _config._nb_partitions  = ( ( volume_per_pass* _config._nb_partitions_in_parallel) / _config._max_memory ) + 1;

        //printf("nb passes  %i  (nb part %i / %zu)\n",_nb_passes,_nb_partitions,max_open_files);
        //_nb_partitions = max_open_files; break;

        if (_config._nb_partitions >= max_open_files && _config._nb_partitions_in_parallel >1)     { _config._nb_partitions_in_parallel  = _config._nb_partitions_in_parallel /2;  }
        else if (_config._nb_partitions >= max_open_files && _config._nb_partitions_in_parallel == 1)   { _config._nb_passes++;  }
        else                                                                            { break;         }

        //printf("update nb passes  %i  (nb part %i / %zu)\n",_nb_passes,_nb_partitions,max_open_files);
    } while (1);

    //if (_config._nb_partitions < 50 &&  (max_open_files - _config._nb_partitions  > 30) ) _config._nb_partitions += 30; //a hack to have more partitions than 30

    //round nb parti to upper multiple of _nb_partitions_in_parallel if possible
    int  incpart = _config._nb_partitions_in_parallel - _config._nb_partitions % _config._nb_partitions_in_parallel;
    incpart = incpart % _config._nb_partitions_in_parallel;
    if(((int)max_open_files - (int)_config._nb_partitions  > incpart)) _config._nb_partitions+= incpart ;

    //_nb_partitions_in_parallel = 1 ;

    //then put _nbCores_per_partition

    _config._nbCores_per_partition =  _config._nbCores / _config._nb_partitions_in_parallel ; //how to best distrib available cores ?

    // with this formula we'll sometimes use less than _nbCores (maybe enforce _nb_partitions_in_parallel is power of two ?)
    DEBUG (("ConfigurationAlgorithm<span>::execute  _nbCores %zu  _nb_partitions_in_parallel %zu   _nbCores_per_partition %zu  nb part %u nb passes %i\n",
        _config._nbCores, _config._nb_partitions_in_parallel, _config._nbCores_per_partition, _config._nb_partitions, _config._nb_passes
    ));

    assert(_config._nbCores_per_partition > 0);

    /* optimize the number of cached items per partition per core */
    /* add more items to partition cache as long as the total memory of cached items 
     * don't use more than a 10th of the requested memory (accepatble overhead) 
     * detail of the memory usage of items in partitions cache:
     * number of cached items * number of partition * number of cores * size of a kmer
     * e.g. (4096 kmers * 8 bytes = 32 KB per partion) * 6000 partitions * 32 cores = 6 GB of buffer */

    uint64_t memoryUsageCachedItems;
    _config._nb_cached_items_per_core_per_part = 1 << 8; // cache at least 256 items (128 here, then * 2 in the next while loop)

    do
    {
        _config._nb_cached_items_per_core_per_part *= 2;
        memoryUsageCachedItems = 1LL * _config._nb_cached_items_per_core_per_part *_config._nb_partitions * _config._nbCores * sizeof(Type); 
    }
    while (memoryUsageCachedItems < _config._max_memory * MBYTE / 10);
        
    DEBUG (("ConfigurationAlgorithm<span>::execute  _config._nb_cached_items_per_core_per_part : %zu ; total memory usage of cached items : %lld MB \n",
        _config._nb_cached_items_per_core_per_part, memoryUsageCachedItems / MBYTE
    ));


    /** We set the config as set. */
    _config._isComputed = true;

    /** We collect some statistics. */
    getInfo()->add (1, _config.getProperties());
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
vector<CountRange> ConfigurationAlgorithm<span>::getSolidityThresholds (IProperties* params)
{
    vector<CountRange> thresholds;

    /** We get the abundance max value. */
    u_int64_t abundanceMax = params->getInt(STR_KMER_ABUNDANCE_MAX);

    /** We split the abundance min string. */
    TokenizerIterator it (params->getStr(STR_KMER_ABUNDANCE_MIN).c_str(), ",");
    for (it.first(); !it.isDone(); it.next())
    {
        CountNumber val = 0;

        if (strcmp(it.item(),"auto")==0)  { val = -1; }
        else                              { val = atoll (it.item());  }

        thresholds.push_back (CountRange(val, abundanceMax));
    }

    return thresholds;
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
	std::vector<bool> ConfigurationAlgorithm<span>::getSolidityCustomVector (IProperties* params)
	{
		std::vector<bool> solidVec;

		std::string solidList = params->getStr (STR_SOLIDITY_CUSTOM);
		for ( std::string::iterator it=solidList.begin(); it!=solidList.end(); ++it)
		{
			solidVec.push_back(*it =='1');
		}
		return solidVec;
	}
	
	
/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

