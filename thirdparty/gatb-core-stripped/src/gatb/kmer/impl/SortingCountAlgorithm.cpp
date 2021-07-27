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

/********************************************************************************/
// We include required definitions
/********************************************************************************/

#include <gatb/kmer/impl/SortingCountAlgorithm.hpp>
#include <gatb/kmer/impl/CountProcessor.hpp>
#include <gatb/kmer/impl/Sequence2SuperKmer.hpp>
#include <gatb/kmer/impl/ConfigurationAlgorithm.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>
#include <gatb/kmer/impl/PartitionsCommand.hpp>
#include <gatb/kmer/impl/RepartitionAlgorithm.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/tools/collections/impl/IterableHelpers.hpp>
#include <cmath>

#define DEBUG(a)  //printf a

/********************************************************************************/
// We use the required packages
/********************************************************************************/
using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::bank;
using namespace gatb::core::bank::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::storage::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::tools::math;

using namespace gatb::core::kmer::impl;

/********************************************************************************/

//#define PROTO_COMP

#ifdef PROTO_COMP
    #define PartitionCacheType  PartitionCacheSorted
    #define STORAGE_TYPE  STORAGE_COMPRESSED_FILE
#else
    #define PartitionCacheType  PartitionCache
    #define STORAGE_TYPE  STORAGE_FILE
#endif

/********************************************************************************/
namespace gatb  {  namespace core  {   namespace kmer  {   namespace impl {
/********************************************************************************/

/********************************************************************************/
static const char* progressFormat0 = "DSK: counting kmers                    ";
static const char* progressFormat1 = "DSK: Pass %d/%d, Step 1: partitioning    ";
static const char* progressFormat2 = "DSK: Pass %d/%d, Step 2: counting kmers  ";
static const char* progressFormat4 = "DSK: nb solid kmers found : %-9ld  ";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
SortingCountAlgorithm<span>::SortingCountAlgorithm (IProperties* params)
  : Algorithm("dsk", -1, params),
    _bank(0), _repartitor(0),
    _progress (0), _tmpPartitionsStorage(0), _tmpPartitions(0), _storage(0),_superKstorage(0)
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
SortingCountAlgorithm<span>::SortingCountAlgorithm (IBank* bank, IProperties* params)
  : Algorithm("dsk", -1, params),
    _bank(0), _repartitor(0),
    _progress (0),_tmpPartitionsStorage(0), _tmpPartitions(0), _storage(0),_superKstorage(0)
{
    setBank (bank);
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
SortingCountAlgorithm<span>::SortingCountAlgorithm (
    IBank*                  bank,
    const Configuration&    config,
    Repartitor*             repartitor,
    vector<CountProcessor*> processors,
	tools::misc::IProperties* params

)
  : Algorithm("dsk", config._nbCores, params),
    _config(config), _bank(0), _repartitor(0),
    _progress (0),_tmpPartitionsStorage(0), _tmpPartitions(0), _storage(0),_superKstorage(0)
{
    setBank       (bank);
    setRepartitor (repartitor);

    for (size_t i=0; i<processors.size(); i++)  {  addProcessor  (processors[i]); }
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
SortingCountAlgorithm<span>::~SortingCountAlgorithm ()
{
    setBank                 (0);
    setRepartitor           (0);
    setProgress             (0);
 //   setPartitionsStorage    (0);
 //   setPartitions           (0);
    setStorage              (0);

    for (size_t i=0; i<_processors.size(); i++)  { _processors[i]->forget(); }
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
SortingCountAlgorithm<span>& SortingCountAlgorithm<span>::operator= (const SortingCountAlgorithm& s)
{
    if (this != &s)
    {
        _config = s._config;

        setBank                 (s._bank);
        setRepartitor           (s._repartitor);
        setProgress             (s._progress);
        setPartitionsStorage    (s._tmpPartitionsStorage);
	    setPartitions           (s._tmpPartitions);
		_superKstorage = s._superKstorage;
        setStorage              (s._storage);
    }
    return *this;
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
IOptionsParser* SortingCountAlgorithm<span>::getOptionsParser (bool mandatory)
{
    IOptionsParser* parser = new OptionsParser ("kmer count");

    string abundanceMax = Stringify::format("%ld", std::numeric_limits<CountNumber>::max());

    parser->push_back (new OptionOneParam (STR_URI_INPUT,         "reads file", mandatory ));
    parser->push_back (new OptionOneParam (STR_KMER_SIZE,         "size of a kmer",                                 false, "31"    ));
    parser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MIN,"min abundance threshold for solid kmers",        false, "2"     ));
    parser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MAX,"max abundance threshold for solid kmers",        false, abundanceMax));
    parser->push_back (new OptionOneParam (STR_KMER_ABUNDANCE_MIN_THRESHOLD,"min abundance hard threshold (only used when min abundance is \"auto\")",false, "2"));
    parser->push_back (new OptionOneParam (STR_HISTOGRAM_MAX,     "max number of values in kmers histogram",        false, "10000"));
    parser->push_back (new OptionOneParam (STR_SOLIDITY_KIND,     "way to compute counts of several files (sum, min, max, one, all, custom)",false, "sum"));
	parser->push_back (new OptionOneParam (STR_SOLIDITY_CUSTOM,   "when solidity-kind is custom, specifies list of files where kmer must be present",false, ""));
    parser->push_back (new OptionOneParam (STR_MAX_MEMORY,        "max memory (in MBytes)",                         false, "5000"));
    parser->push_back (new OptionOneParam (STR_MAX_DISK,          "max disk   (in MBytes)",                         false, "0"));
    parser->push_back (new OptionOneParam (STR_URI_SOLID_KMERS,   "output file for solid kmers (only when constructing a graph)", false));
    parser->push_back (new OptionOneParam (STR_URI_OUTPUT,        "output file",                                    false));
    parser->push_back (new OptionOneParam (STR_URI_OUTPUT_DIR,    "output directory",                               false, "."));
    parser->push_back (new OptionOneParam (STR_URI_OUTPUT_TMP,    "output directory for temporary files",           false, "."));
    parser->push_back (new OptionOneParam (STR_STORAGE_TYPE,      "storage type of kmer counts ('file')", false, "file"  ));
	parser->push_back (new OptionOneParam (STR_HISTO2D,"compute the 2D histogram (with first file = genome, remaining files = reads)",false,"0"));
	parser->push_back (new OptionOneParam (STR_HISTO,"output the kmer abundance histogram",false,"0"));

    IOptionsParser* devParser = new OptionsParser ("kmer count, advanced performance tweaks");

    devParser->push_back (new OptionOneParam (STR_MINIMIZER_TYPE,    "minimizer type (0=lexi, 1=freq)",                false, "0"));
    devParser->push_back (new OptionOneParam (STR_MINIMIZER_SIZE,    "size of a minimizer",                            false, "10"));
    devParser->push_back (new OptionOneParam (STR_REPARTITION_TYPE,  "minimizer repartition (0=unordered, 1=ordered)", false, "0"));
    parser->push_back (devParser);

    return parser;
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
IProperties* SortingCountAlgorithm<span>::getDefaultProperties ()
{
    IOptionsParser* parser = getOptionsParser (true);
    LOCAL (parser);
    return parser->getDefaultProperties();
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
ICountProcessor<span>* SortingCountAlgorithm<span>::getDefaultProcessor (
    tools::misc::IProperties*       params,
    tools::storage::impl::Storage*  dskStorage,
    tools::storage::impl::Storage*  otherStorage
)
{
	
	std::string histo2Dstorage_filename;
	// build histo2D filename
	int using_histo_2D = params->get(STR_HISTO2D) ? params->getInt(STR_HISTO2D) : 0 ;
	if(using_histo_2D)
	{
		if(params->get(STR_URI_OUTPUT))
		{
			std::string uri_input = params->getStr(STR_URI_OUTPUT);
			histo2Dstorage_filename =  uri_input + ".histo2D";
		}
		else
		if(params->get(STR_URI_INPUT))
		{
			std::string uri_input = params->getStr(STR_URI_INPUT);
			std::string delimiter = ",";
			std::string firstbankname = uri_input.substr(0, uri_input.find(delimiter));
			histo2Dstorage_filename =  system::impl::System::file().getBaseName(firstbankname) + ".histo2D";
		}
		else if(params->get(STR_URI_FILE))
		{
			std::string uri_input = params->getStr(STR_URI_FILE);
			std::string delimiter = ",";
			std::string firstbankname = uri_input.substr(0, uri_input.find(delimiter));
			histo2Dstorage_filename =  system::impl::System::file().getBaseName(firstbankname) + ".histo2D";
		}
		else
		{
			histo2Dstorage_filename = "histo2D_resultfile";
		}

	}
	//histo1D filename
	std::string histo1Dstorage_filename;

	int using_histo_1D = params->get(STR_HISTO) ? params->getInt(STR_HISTO) : 0 ;
	if(using_histo_1D)
	{
		if(params->get(STR_URI_OUTPUT))
		{
			std::string uri_input = params->getStr(STR_URI_OUTPUT);
			histo1Dstorage_filename =  uri_input + ".histo";
		}
		else
			if(params->get(STR_URI_INPUT))
			{
				std::string uri_input = params->getStr(STR_URI_INPUT);
				std::string delimiter = ",";
				std::string firstbankname = uri_input.substr(0, uri_input.find(delimiter));
				histo1Dstorage_filename =  system::impl::System::file().getBaseName(firstbankname) + ".histo";
			}
			else if(params->get(STR_URI_FILE))
			{
				std::string uri_input = params->getStr(STR_URI_FILE);
				std::string delimiter = ",";
				std::string firstbankname = uri_input.substr(0, uri_input.find(delimiter));
				histo1Dstorage_filename =  system::impl::System::file().getBaseName(firstbankname) + ".histo";
			}
			else
			{
				histo1Dstorage_filename = "histo_resultfile";
			}
		
	}
	
	
    CountProcessor* result = 0;

    if (otherStorage == 0)  { otherStorage = dskStorage; }

    if (params==0 || dskStorage==0 || otherStorage==0)  { throw Exception ("Bad parameters in SortingCountAlgorithm<span>::getDefaultProcessor"); }

    /** The default count processor is defined as the following chain :
     *      1) histogram
     *      2) solidity filter
     *      3) if solidity filter passed, dump to file system
     */
    result = new CountProcessorChain<span> (

        new CountProcessorHistogram<span> (
            & otherStorage->getGroup("histogram"),
            params->getInt(STR_HISTOGRAM_MAX),
            params->getInt(STR_KMER_ABUNDANCE_MIN_THRESHOLD),
			using_histo_2D,
			using_histo_1D,
			histo2Dstorage_filename,
			histo1Dstorage_filename
        ),

        CountProcessorSolidityFactory<span>::create (*params),

        new CountProcessorDump     <span> (
            dskStorage->getGroup("dsk"),
            params->getInt(STR_KMER_SIZE)
        ),
        NULL
    );

    /** We set some name. */
    result->setName ("dsk");

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

/** We need a specific count processor proxy. */
template<size_t span>
class CountProcessorCustomProxy : public CountProcessorProxy<span>
{
    public:
    CountProcessorCustomProxy (ICountProcessor<span>* cutoffProcessor, ICountProcessor<span>* dskProcessor)
        : CountProcessorProxy<span>(cutoffProcessor), _cutoffProcessor(cutoffProcessor), _dskProcessor(dskProcessor) {}

    /** \copydoc ICountProcessor<span>::end */
    void endPass (size_t passId)
    {
        /** We call the parent method. */
        CountProcessorProxy<span>::endPass (passId);

        /** Now, we have the cutoffs information, and we can put it as abundance min of the dsk processor. */
        if (CountProcessorCutoff<span>* cutoffProc = dynamic_cast<CountProcessorCutoff<span>*> (_cutoffProcessor))
        {
            if (CountProcessorSolidityInfo* info = _dskProcessor->template get<CountProcessorSolidityInfo> ())
            {
                info->setAbundanceMin (cutoffProc->getCutoffs());
            }
        }
    }

    private:
        ICountProcessor<span>* _cutoffProcessor;
        ICountProcessor<span>* _dskProcessor;
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
vector<ICountProcessor<span>*> SortingCountAlgorithm<span>::getDefaultProcessorVector (
    Configuration&  config,
    IProperties*    params,
    Storage*        dskStorage,
    Storage*        otherStorage
)
{

    vector<ICountProcessor<span>*> result;

    ICountProcessor<span>* dskProcessor = getDefaultProcessor (params, dskStorage, otherStorage);

    /** Now, we define the vector of count processors to be given to the SortingCountAlgorithm.
     * The choice depends on the presence of "auto" min abundance in the configuration. */
    bool foundAuto = false;
    for (size_t i=0; !foundAuto && i<config._abundance.size(); i++) {  foundAuto = config._abundance[i].getBegin() == -1;  }

    if (foundAuto)
    {
        /** We create a cutoff count processor to compute cutoffs. */
        ICountProcessor<span>* cutoffProcessor = 0;

        /** We look what kind of solidity is required. If it is in [sum,min,max], we consider that we have only one bank
         * to compute the histogram for, even if N banks are provided. In such a case, the histogram is computed for the
         * concatenation of the N banks. */
        switch (config._solidityKind)
        {
        case KMER_SOLIDITY_MIN:
        case KMER_SOLIDITY_MAX:
        case KMER_SOLIDITY_SUM:
            cutoffProcessor = new CountProcessorCutoff<span> (1);
            break;

        case KMER_SOLIDITY_ONE:
		case KMER_SOLIDITY_CUSTOM:
        case KMER_SOLIDITY_ALL:
            /** We create a cutoff count processor to compute cutoffs of the N banks. */
            cutoffProcessor = new CountProcessorCutoff<span> (config._nb_banks);
            break;

        default:
            break;
        }

        if (cutoffProcessor == 0)  { throw Exception ("Unable to configure count processor due to bad solidity kind %d", config._solidityKind); }

        /** We encapsulate both cutoff and dsk processors in a single one that link them. */
        ICountProcessor<span>* proxyCutoff = new CountProcessorCustomProxy<span> (cutoffProcessor, dskProcessor);
        proxyCutoff->setName("cutoffs_auto");

        result.push_back (proxyCutoff);
        result.push_back (dskProcessor);
    }
    else
    {
        result.push_back (dskProcessor);
    }

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
template<size_t span>
void SortingCountAlgorithm<span>::configure ()
{
    DEBUG (("SortingCountAlgorithm<span>::configure  BEGIN  _bank=%p  _config.isComputed=%d  _repartitor=%p  \n",
        _bank, _config._isComputed, _repartitor
    ));

    /** We check that the bank is ok, otherwise we build one. */
    if (_bank == 0)    {  setBank (Bank::open (getInput()->getStr(STR_URI_INPUT)));  }

    /** We may have to create a default storage. */
    Storage* storage = 0;
    if (_repartitor==0 || _processors.size() == 0)
    {
        string output = getInput()->get(STR_URI_OUTPUT) ?
            getInput()->getStr(STR_URI_OUTPUT)   :
            (getInput()->getStr(STR_URI_OUTPUT_DIR) + "/" + system::impl::System::file().getBaseName (_bank->getId()));

        /* create output dir if it doesn't exist */
        if(!System::file().doesExist(getInput()->getStr(STR_URI_OUTPUT_DIR))){
            int ok = System::file().mkdir(getInput()->getStr(STR_URI_OUTPUT_DIR), 0755);
            if(ok != 0){
                throw Exception ("Error: can't create output directory");
            }
        }

        string storage_type = getInput()->getStr(STR_STORAGE_TYPE);
        if (storage_type == "file")
            _storage_type = tools::storage::impl::STORAGE_FILE;
        else
        {std::cout << "Error: unknown storage type specified: " << storage_type << std::endl; exit(1); }

        storage = StorageFactory(_storage_type).create (output, true, false); //// this is the storage for the output (kmer counts)
    }

    /** In case the storage is created in this method, we need to keep an eye on it. */
    setStorage (storage);

    /** We check that the configuration is ok, otherwise we build one. */
    if (_config._isComputed == false)
    {
        ConfigurationAlgorithm<span> configAlgo (_bank, getInput());
        configAlgo.execute();
        _config = configAlgo.getConfiguration();
 
        /* remember configuration details (e.g. number of passes, partitions). useful for bcalm. */
        storage->getGroup(configAlgo.getName()).setProperty("xml", string("\n") + configAlgo.getInfo()->getXML());
   }

    /** We check that the minimizers hash function is ok, otherwise we build one. */
    if (_repartitor == 0)
    {
        RepartitorAlgorithm<span> repart (
                _bank, 
                storage->getGroup("minimizers"), 
                _config,
                getInput()->get(STR_NB_CORES) ? getInput()->getInt(STR_NB_CORES) : 0
                );
        repart.execute ();
        setRepartitor (new Repartitor(storage->getGroup("minimizers")));
    }

	
	//out file name for histo2D
	string output_histo2Dname  = getInput()->get(STR_URI_OUTPUT) ?
	getInput()->getStr(STR_URI_OUTPUT) + ".histo2D"   :
	(getInput()->getStr(STR_URI_OUTPUT_DIR) + "/" + system::impl::System::file().getBaseName (_bank->getIdNb(0))) + ".histo2D";
	
	//printf("out histo name %s  \n",output_histo2Dname.c_str());

	/** We check that the processor is ok, otherwise we build one. */
    if (_processors.size() == 0)  { _processors = getDefaultProcessorVector(_config, getInput(), storage, storage);  };

	
	/** We check if options are compatible with histo2D */
	int using_histo_2D = getInput()->get(STR_HISTO2D) ? getInput()->getInt(STR_HISTO2D) : 0 ;
	if(using_histo_2D)
	{

		getInput()->setStr(STR_SOLIDITY_KIND, "all");
		_config._solidityKind = KMER_SOLIDITY_ALL;
	
		//check number of banks, must be greater than 1
		int nbanks = _bank->getBanks().size();

		if( nbanks <2 )
		{
			fprintf(stderr,"There must be at least 2 input banks when using -histo2D \n");
			exit(1);
		}
		
	}
	
    DEBUG (("SortingCountAlgorithm<span>::configure  END  _bank=%p  _config.isComputed=%d  _repartitor=%p  storage=%p\n",
        _bank, _config._isComputed, _repartitor, storage
    ));
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
void SortingCountAlgorithm<span>::execute ()
{
    /*************************************************************/
    /*                       CONFIGURATION                       */
    /*************************************************************/

    /** We configure all required objects (bank, configuration, repartitor, count processor). */
    configure ();

    /** We create the sequences iterator. */
    Iterator<Sequence>* itSeq = _bank->iterator();
    LOCAL (itSeq);

    /** We configure the progress bar. Note that we create a ProgressSynchro since this progress bar
     * may me modified by several threads at the same time. */
    size_t nbIterations = (1 + _processors.size()) * _config._volume * MBYTE / sizeof(Type);
    setProgress (new ProgressSynchro (
        createIteratorListener (nbIterations, progressFormat0),
        System::thread().newSynchronizer())
    );
    _progress->init ();

#ifdef NONCANONICAL
    std::cerr << std::endl << "NOTICE: This version was compiled to perform non-canonical kmer counting." << std::endl;
#endif

    if (_config._kmerSize <= 2)
    {
        std::cout << "k-mer counting with k<=2 is not supported" << std::endl; // it's buggy. try it with https://github.com/GATB/dsk/blob/master/test/shortread.fasta, you will see only 5 kmers returned
        exit(1);
    }

    /** We create the PartiInfo instance. */
    PartiInfo<5> pInfo (_config._nb_partitions, _config._minim_size);

    /** We notify the count processor about the start of the main loop. */
    for (size_t i=0; i<_processors.size(); i++)  {  _processors[i]->begin (_config); }

    /*************************************************************/
    /*                         MAIN LOOP                         */
    /*************************************************************/
    /** We loop N times the bank. For each pass, we will consider a subset of the whole kmers set of the bank. */
    for (size_t current_pass=0; current_pass < _config._nb_passes; current_pass++)
    {
        DEBUG (("SortingCountAlgorithm<span>::execute  pass [%ld,%d] \n", current_pass+1, _config._nb_passes));

        pInfo.clear();

        /** 1) We fill the partition files. */
        fillPartitions (current_pass, itSeq, pInfo);

        /** 2) We fill the kmers solid file from the partition files. */
        fillSolidKmers (current_pass, pInfo);
    }

    /** We notify the count processor about the stop of the main loop. */
    for (size_t i=0; i<_processors.size(); i++)  {  _processors[i]->end (); }

    /** We update the progress information. */
    for (size_t i=0; i<_processors.size(); i++)
    {
        CountProcessorDump<span>* processorDump = _processors[i]->template get <CountProcessorDump<span> > ();
        if (processorDump != 0) {  _progress->setMessage (Stringify::format(progressFormat4, processorDump->getNbItems())); }
    }

    _progress->finish ();

	
//	pInfo.printInfo();
	

    /** We want to remove physically the partitions. */
	if(_config._solidityKind != KMER_SOLIDITY_SUM)
     _tmpPartitions->remove ();

	u_int64_t totaltmp, biggesttmp, smallesttmp;
	float meantmp;
	if(_config._solidityKind == KMER_SOLIDITY_SUM)
		_superKstorage->getFilesStats(totaltmp,biggesttmp,smallesttmp, meantmp);


	if(_superKstorage!=0)
	{
		delete _superKstorage; //delete files and containing dir
		_superKstorage =0;
	}
	
    /*************************************************************/
    /*                         STATISTICS                        */
    /*************************************************************/

    /** We gather some statistics. */
    if (_bankStats.sequencesNb > 0)
    {
        getInfo()->add (1, "bank");
        getInfo()->add (2, "bank_uri",          "%s",   _bank->getId().c_str());
        getInfo()->add (2, "bank_size",         "%lld", _bank->getSize());
        getInfo()->add (2, "bank_total_nt",     "%lld", _bankStats.sequencesTotalLength);
        getInfo()->add (2, "sequences");
        getInfo()->add (3, "seq_number",        "%ld",  _bankStats.sequencesNb);
        getInfo()->add (3, "seq_size_min",      "%ld",  _bankStats.sequencesMinLength);
        getInfo()->add (3, "seq_size_max",      "%ld",  _bankStats.sequencesMaxLength);
        getInfo()->add (3, "seq_size_mean",     "%.1f", _bankStats.getSeqMean());
        getInfo()->add (3, "seq_size_deviation","%.1f", _bankStats.getSeqDeviation());
        getInfo()->add (2, "kmers");
        getInfo()->add (3, "kmers_nb_valid",   "%lld", _bankStats.kmersNbValid);
        getInfo()->add (3, "kmers_nb_invalid", "%lld", _bankStats.kmersNbInvalid);
		
		
    }


	
	u_int64_t nbtotalsuperk = pInfo.getNbSuperKmerTotal();
	u_int64_t nbtotalk = pInfo.getNbKmerTotal();
	
    getInfo()->add (1, "stats");
	
	getInfo()->add (2, "temp_files");
	getInfo()->add (3, "nb_superkmers","%lld",nbtotalsuperk);
	getInfo()->add (3, "avg_superk_length","%.2f",(nbtotalk/(float) nbtotalsuperk));
	getInfo()->add (3, "minimizer_density","%.2f",(nbtotalsuperk/(float)nbtotalk)*(_config._kmerSize - _config._minim_size +2));
	
	if(_config._solidityKind == KMER_SOLIDITY_SUM)
	{
		getInfo()->add (3, "total_size_(MB)","%lld",totaltmp/1024LL/1024LL);
		getInfo()->add (3, "tmp_file_biggest_(MB)","%lld",biggesttmp/1024LL/1024LL);
		getInfo()->add (3, "tmp_file_smallest_(MB)","%lld",smallesttmp/1024LL/1024LL);
		getInfo()->add (3, "tmp_file_mean_(MB)","%.1f",meantmp/1024LL/1024LL);
	}
    /** We dump information about count processors. */
    if (_processors.size()==1)  {  getInfo()->add (2, _processors[0]->getProperties()); }
    else
    {
        for (size_t i=0; i<_processors.size(); i++)
        {
            getInfo()->add (2, _processors[i]->getName());
            getInfo()->add (3, _processors[i]->getProperties());
        }
    }

    _fillTimeInfo /= getDispatcher()->getExecutionUnitsNumber();
    getInfo()->add (2, _fillTimeInfo.getProperties("fillsolid_time"));

    getInfo()->add (1, getTimeInfo().getProperties("time"));
}

/********************************************************************************/
/* This functor class takes a Sequence as input, splits it into super kmers and
 * serialize them into partitions.
 *
 * A superkmer is dumped into a partition 'p' chosen by some hash code of the minimizer
 * of the superkmer. Such a hash code can be computed in several way; actually, we use
 * a lookup table that has computed the minimizers distribution on a subset of the
 * processed bank.
 */
	
template<size_t span, bool newSuperKmerStorage=true> //
class FillPartitions : public Sequence2SuperKmer<span>
{
public:
	/** Shortcut. */
	typedef typename Sequence2SuperKmer<span>::Type            Type;
	typedef typename Sequence2SuperKmer<span>::Model           Model;
	typedef typename Model::Kmer                          KmerType;
	typedef typename Kmer<span>::SuperKmer                SuperKmer;
	
	/** */
	void processSuperkmer (SuperKmer& superKmer)
	{
		if ((superKmer.minimizer % this->_nbPass) == this->_pass && superKmer.isValid()) //check if falls into pass
		{
			/** We get the hash code for the current miminizer.
			 * => this will give us the partition where to dump the superkmer. */
			size_t p = this->_repartition (superKmer.minimizer);
			
			/** We save the superkmer into the right partition. */
			if(newSuperKmerStorage)
				superKmer.save (_superkmerFiles,p);
			else
				superKmer.save ((this->_partition)[p]);
			
			
			//for debug purposes
			_local_pInfo.incSuperKmer_per_minimBin (superKmer.minimizer, superKmer.size()); //tocheck
			
			/*********************************************/
			/** Now, we compute statistics about kxmers. */
			/*********************************************/
			
			Type radix_kxmer_forward ,radix_kxmer ;
			bool prev_which = superKmer[0].which();
			size_t kx_size =0;
			
			radix_kxmer_forward = getHeavyWeight (superKmer[0].value());
			
			for (size_t ii=1 ; ii < superKmer.size(); ii++)
			{
				//compute here stats on  kx mer
				//tant que tai <= xmer et which kmer[ii] == which kmer [ii-1] --> cest un kxmer
				//do the same in sampling : gives ram estimation
				if (superKmer[ii].which() != prev_which || kx_size >= this->_kx) // kxmer_size = 1 //cost should diminish with larger kxmer
				{
					//output kxmer size kx_size,radix_kxmer
					//kx mer is composed of _superKp[ii-1] _superKp[ii-2] .. _superKp[ii-n] with nb elems  n  == kxmer_size +1  (un seul kmer ==k+0)
					if(prev_which)
					{
						radix_kxmer = radix_kxmer_forward;
					}
					else // si revcomp, le radix du kxmer est le debut du dernier kmer
					{
						radix_kxmer = getHeavyWeight (superKmer[ii-1].value());
					}
					
					this->_local_pInfo.incKmer_and_rad (p, radix_kxmer.getVal(), kx_size); //nb of superkmer per x per parti per radix
					
					radix_kxmer_forward =  getHeavyWeight (superKmer[ii].value());
					kx_size =0;
				}
				else
				{
					kx_size++;
				}
				
				prev_which = superKmer[ii].which() ;
			}
			
			//record last kx mer
			if(prev_which)
			{
				radix_kxmer = radix_kxmer_forward;
			}
			else // si revcomp, le radix du kxmer est le debut du dernier kmer
			{
				radix_kxmer =  getHeavyWeight (superKmer[superKmer.size()-1].value());
			}
			
			this->_local_pInfo.incKmer_and_rad(p, radix_kxmer.getVal(),kx_size );
			
			/** We update progression information. */
			this->_nbWrittenKmers += superKmer.size();
		}
	}
	
	/** Constructor. */
	FillPartitions (
					Model&             model,
					size_t             nbPasses,
					size_t             currentPass,
					size_t             nbPartitions,
					size_t             nbCacheItems,
					IteratorListener*  progress,
					BankStats&         bankStats,
					Partition<Type>*   partition,
					Repartitor&        repartition,
					PartiInfo<5>&      pInfo,
					SuperKmerBinFiles* superKstorage
					)
	:   Sequence2SuperKmer<span> (model, nbPasses, currentPass, nbPartitions, progress, bankStats),
	_kx(4),
	_extern_pInfo(pInfo) , _local_pInfo(nbPartitions,model.getMmersModel().getKmerSize()),
	_repartition (repartition)
	,_partition (*partition, nbCacheItems, 0), _superkmerFiles(superKstorage, nbCacheItems* sizeof(Type))
	{
		_mask_radix.setVal((int64_t) 255);
		_mask_radix = _mask_radix << ((this->_kmersize - 4)*2); //get first 4 nt  of the kmers (heavy weight)
	}
	
	
	/** Destructor. */
	virtual ~FillPartitions ()
	{
		//printf("destruc fillparti _superkmerFiles %p \n",_superkmerFiles);
		
		//add to global parti_info
		_extern_pInfo.add_sync(_local_pInfo);
	}
	
private:
	
	size_t        _kx;
	PartiInfo<5>& _extern_pInfo;
	PartiInfo<5>  _local_pInfo;
	Type          _mask_radix;
	Repartitor&   _repartition;
	
	/** Shared resources (must support concurrent accesses). */
	PartitionCacheType <Type> _partition;
		CacheSuperKmerBinFiles _superkmerFiles;
	
	
	
	Type getHeavyWeight (const Type& kmer) const  {  return (kmer & this->_mask_radix) >> ((this->_kmersize - 4)*2);  }
};
	
	
//specialization when newSuperKmerStorage = false
template<size_t span>
class FillPartitions <span, false>: public Sequence2SuperKmer<span>
{
public:
    /** Shortcut. */
    typedef typename Sequence2SuperKmer<span>::Type            Type;
    typedef typename Sequence2SuperKmer<span>::Model           Model;
    typedef typename Model::Kmer                          KmerType;
    typedef typename Kmer<span>::SuperKmer                SuperKmer;

    /** */
    void processSuperkmer (SuperKmer& superKmer)
    {
        if ((superKmer.minimizer % this->_nbPass) == this->_pass && superKmer.isValid()) //check if falls into pass
        {
            /** We get the hash code for the current miminizer.
             * => this will give us the partition where to dump the superkmer. */
            size_t p = this->_repartition (superKmer.minimizer);

            /** We save the superkmer into the right partition. */
//			if(newSuperKmerStorage)
//				superKmer.save (_superkmerFiles,p);
//			else
				superKmer.save ((this->_partition)[p]);

	
			//for debug purposes
			_local_pInfo.incSuperKmer_per_minimBin (superKmer.minimizer, superKmer.size()); //tocheck

            /*********************************************/
            /** Now, we compute statistics about kxmers. */
            /*********************************************/

            Type radix_kxmer_forward ,radix_kxmer ;
            bool prev_which = superKmer[0].which();
            size_t kx_size =0;

            radix_kxmer_forward = getHeavyWeight (superKmer[0].value());

            for (size_t ii=1 ; ii < superKmer.size(); ii++)
            {
                //compute here stats on  kx mer
                //tant que tai <= xmer et which kmer[ii] == which kmer [ii-1] --> cest un kxmer
                //do the same in sampling : gives ram estimation
                if (superKmer[ii].which() != prev_which || kx_size >= this->_kx) // kxmer_size = 1 //cost should diminish with larger kxmer
                {
                    //output kxmer size kx_size,radix_kxmer
                    //kx mer is composed of _superKp[ii-1] _superKp[ii-2] .. _superKp[ii-n] with nb elems  n  == kxmer_size +1  (un seul kmer ==k+0)
                    if(prev_which)
                    {
                        radix_kxmer = radix_kxmer_forward;
                    }
                    else // si revcomp, le radix du kxmer est le debut du dernier kmer
                    {
                        radix_kxmer = getHeavyWeight (superKmer[ii-1].value());
                    }

                    this->_local_pInfo.incKmer_and_rad (p, radix_kxmer.getVal(), kx_size); //nb of superkmer per x per parti per radix

                    radix_kxmer_forward =  getHeavyWeight (superKmer[ii].value());
                    kx_size =0;
                }
                else
                {
                    kx_size++;
                }

                prev_which = superKmer[ii].which() ;
            }

            //record last kx mer
            if(prev_which)
            {
                radix_kxmer = radix_kxmer_forward;
            }
            else // si revcomp, le radix du kxmer est le debut du dernier kmer
            {
                radix_kxmer =  getHeavyWeight (superKmer[superKmer.size()-1].value());
            }

            this->_local_pInfo.incKmer_and_rad(p, radix_kxmer.getVal(),kx_size );

            /** We update progression information. */
            this->_nbWrittenKmers += superKmer.size();
        }
    }

    /** Constructor. */
    FillPartitions (
        Model&             model,
        size_t             nbPasses,
        size_t             currentPass,
        size_t             nbPartitions,
        size_t             nbCacheItems,
        IteratorListener*  progress,
        BankStats&         bankStats,
        Partition<Type>*   partition,
        Repartitor&        repartition,
        PartiInfo<5>&      pInfo,
		SuperKmerBinFiles* superKstorage
    )
    :   Sequence2SuperKmer<span> (model, nbPasses, currentPass, nbPartitions, progress, bankStats),
        _kx(4),
        _extern_pInfo(pInfo) , _local_pInfo(nbPartitions,model.getMmersModel().getKmerSize()),
        _repartition (repartition)
	    ,_partition (*partition, nbCacheItems, 0)/*, _superkmerFiles(superKstorage,nbCacheItems* sizeof(Type))*/
    {
        _mask_radix.setVal((int64_t) 255);
        _mask_radix = _mask_radix << ((this->_kmersize - 4)*2); //get first 4 nt  of the kmers (heavy weight)
    }


    /** Destructor. */
    virtual ~FillPartitions ()
    {
		//printf("destruc fillparti _superkmerFiles %p \n",_superkmerFiles);

        //add to global parti_info
        _extern_pInfo.add_sync(_local_pInfo);
    }

private:

    size_t        _kx;
    PartiInfo<5>& _extern_pInfo;
    PartiInfo<5>  _local_pInfo;
    Type          _mask_radix;
    Repartitor&   _repartition;
	
    /** Shared resources (must support concurrent accesses). */
    PartitionCacheType <Type> _partition;
//	CacheSuperKmerBinFiles _superkmerFiles;
	
	

    Type getHeavyWeight (const Type& kmer) const  {  return (kmer & this->_mask_radix) >> ((this->_kmersize - 4)*2);  }
};

	
	//specialization when newSuperKmerStorage = true
	template<size_t span>
	class FillPartitions<span, true> : public Sequence2SuperKmer<span>
	{
	public:
		/** Shortcut. */
		typedef typename Sequence2SuperKmer<span>::Type            Type;
		typedef typename Sequence2SuperKmer<span>::Model           Model;
		typedef typename Model::Kmer                          KmerType;
		typedef typename Kmer<span>::SuperKmer                SuperKmer;
		
		/** */
		void processSuperkmer (SuperKmer& superKmer)
		{
			if ((superKmer.minimizer % this->_nbPass) == this->_pass && superKmer.isValid()) //check if falls into pass
			{
				/** We get the hash code for the current miminizer.
				 * => this will give us the partition where to dump the superkmer. */
				size_t p = this->_repartition (superKmer.minimizer);
				
				/** We save the superkmer into the right partition. */
				superKmer.save (_superkmerFiles,p);

				//for debug purposes
				_local_pInfo.incSuperKmer_per_minimBin (superKmer.minimizer, superKmer.size()); //tocheck
				
				/*********************************************/
				/** Now, we compute statistics about kxmers. */
				/*********************************************/
				
				Type radix_kxmer_forward ,radix_kxmer ;
				bool prev_which = superKmer[0].which();
				size_t kx_size =0;
				
				radix_kxmer_forward = getHeavyWeight (superKmer[0].value());
				
				for (size_t ii=1 ; ii < superKmer.size(); ii++)
				{
					//compute here stats on  kx mer
					//tant que tai <= xmer et which kmer[ii] == which kmer [ii-1] --> cest un kxmer
					//do the same in sampling : gives ram estimation
					if (superKmer[ii].which() != prev_which || kx_size >= this->_kx) // kxmer_size = 1 //cost should diminish with larger kxmer
					{
						//output kxmer size kx_size,radix_kxmer
						//kx mer is composed of _superKp[ii-1] _superKp[ii-2] .. _superKp[ii-n] with nb elems  n  == kxmer_size +1  (un seul kmer ==k+0)
						if(prev_which)
						{
							radix_kxmer = radix_kxmer_forward;
						}
						else // si revcomp, le radix du kxmer est le debut du dernier kmer
						{
							radix_kxmer = getHeavyWeight (superKmer[ii-1].value());
						}
						
						this->_local_pInfo.incKmer_and_rad (p, radix_kxmer.getVal(), kx_size); //nb of superkmer per x per parti per radix
						
						radix_kxmer_forward =  getHeavyWeight (superKmer[ii].value());
						kx_size =0;
					}
					else
					{
						kx_size++;
					}
					
					prev_which = superKmer[ii].which() ;
				}
				
				//record last kx mer
				if(prev_which)
				{
					radix_kxmer = radix_kxmer_forward;
				}
				else // si revcomp, le radix du kxmer est le debut du dernier kmer
				{
					radix_kxmer =  getHeavyWeight (superKmer[superKmer.size()-1].value());
				}
				
				this->_local_pInfo.incKmer_and_rad(p, radix_kxmer.getVal(),kx_size );
				
				/** We update progression information. */
				this->_nbWrittenKmers += superKmer.size();
			}
		}
		
		/** Constructor. */
		FillPartitions (
						Model&             model,
						size_t             nbPasses,
						size_t             currentPass,
						size_t             nbPartitions,
						size_t             nbCacheItems,
						IteratorListener*  progress,
						BankStats&         bankStats,
						Partition<Type>*   partition,
						Repartitor&        repartition,
						PartiInfo<5>&      pInfo,
						SuperKmerBinFiles* superKstorage
						)
		:   Sequence2SuperKmer<span> (model, nbPasses, currentPass, nbPartitions, progress, bankStats),
		_kx(4),
		_extern_pInfo(pInfo) , _local_pInfo(nbPartitions,model.getMmersModel().getKmerSize()),
		_repartition (repartition)
		,/*_partition (*partition, nbCacheItems, 0),*/ _superkmerFiles(superKstorage,nbCacheItems* sizeof(Type))
		{
			_mask_radix.setVal((int64_t) 255);
			_mask_radix = _mask_radix << ((this->_kmersize - 4)*2); //get first 4 nt  of the kmers (heavy weight)
		}
		
		
		/** Destructor. */
		virtual ~FillPartitions ()
		{
			//printf("destruc fillparti _superkmerFiles %p \n",_superkmerFiles);
			
			//add to global parti_info
			_extern_pInfo.add_sync(_local_pInfo);
		}
		
	private:
		
		size_t        _kx;
		PartiInfo<5>& _extern_pInfo;
		PartiInfo<5>  _local_pInfo;
		Type          _mask_radix;
		Repartitor&   _repartition;
		
		/** Shared resources (must support concurrent accesses). */
		//PartitionCacheType <Type> _partition;
		CacheSuperKmerBinFiles _superkmerFiles;
		
		
		
		Type getHeavyWeight (const Type& kmer) const  {  return (kmer & this->_mask_radix) >> ((this->_kmersize - 4)*2);  }
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
void SortingCountAlgorithm<span>::fillPartitions (size_t pass, Iterator<Sequence>* itSeq, PartiInfo<5>& pInfo)
	{
		TIME_INFO (getTimeInfo(), "fill_partitions");
		
		DEBUG (("SortingCountAlgorithm<span>::fillPartitions  _kmerSize=%d _minim_size=%d \n", _config._kmerSize, _config._minim_size));
		
		_nbKmersPerPartitionPerBank.clear();
		if(_config._solidityKind != KMER_SOLIDITY_SUM)
		{
			/** We delete the previous partitions storage. */
			if (_tmpPartitionsStorage)  { _tmpPartitionsStorage->remove (); }
			
			/** We build the temporary storage name from the output storage name. */
			string tmpStorageName = getInput()->getStr(STR_URI_OUTPUT_TMP) + "/" + System::file().getTemporaryFilename("dsk_partitions");
			
			/** We create the partition files for the current pass. */
			setPartitionsStorage (StorageFactory(STORAGE_TYPE).create (tmpStorageName, true, false));
			setPartitions        (0); // close the partitions first, otherwise new files are opened before  closing parti from previous pass
			setPartitions        ( & (*_tmpPartitionsStorage)().getPartition<Type> ("parts", _config._nb_partitions));
			
		}
		else
		{
			/** We build the temporary storage name from the output storage name. */
			_tmpStorageName_superK = getInput()->getStr(STR_URI_OUTPUT_TMP) + "/" + System::file().getTemporaryFilename("superK_partitions");
			
			
			if(_superKstorage!=0)
			{
				delete _superKstorage;
				_superKstorage =0;
			}
			
			_superKstorage = new SuperKmerBinFiles(_tmpStorageName_superK,"superKparts", _config._nb_partitions, false) ;
			
		}
		/** We update the message of the progress bar. */
		_progress->setMessage (Stringify::format(progressFormat1, pass+1, _config._nb_passes));
		
		/** We create a kmer model; using the frequency order if we're in that mode */
		uint32_t* freq_order = NULL;
		
		/** We may have to retrieve the minimizers frequencies computed in the RepartitorAlgorithm. */
		if (_config._minimizerType == 1)  {  freq_order = _repartitor->getMinimizerFrequencies ();  }
		
		Model model( _config._kmerSize, _config._minim_size, typename kmer::impl::Kmer<span>::ComparatorMinimizerFrequencyOrLex(), freq_order);
		
		/** We have to reinit the progress instance since it may have been used by SampleRepart before. */
		_progress->init();
		
		if (_config._solidityKind == KMER_SOLIDITY_SUM) {
			/** We launch the iteration of the sequences iterator with the created
			 * functors. */

			size_t groupSize = 1000;
			bool deleteSynchro = true;

			/** We fill the partitions. Each thread will read synchronously and will
			 * call FillPartitions in a synchronous way (in order to have global
			 * BanksStats correctly computed). */
			getDispatcher()->iterate(
				itSeq,
				FillPartitions<span, true>(
				    model, _config._nb_passes, pass, _config._nb_partitions,
				    _config._nb_cached_items_per_core_per_part, _progress, _bankStats,
				    _tmpPartitions, *_repartitor, pInfo, _superKstorage),
				groupSize, deleteSynchro);

			// GR: close the input bank here with call to finalize
			itSeq->finalize();

			_superKstorage->flushFiles();
			_superKstorage->closeFiles();
		} else {
			/** We may have several input banks instead of a single one. */
			std::vector<Iterator<Sequence>*> itBanks =  itSeq->getComposition();

			/** We first reset the vector holding the kmers number for each partition and for each bank.
			 * It can be seen as the following matrix:
			 *
			 *           part0  part1  part2 ... partJ
			 *   bank0    xxx    xxx    xxx       xxx
			 *   bank1    xxx    xxx    xxx       xxx
			 *    ...
			 *   bankI    xxx    xxx    xxx       xxx
			 *
			 *   Here xxx is the number of items found for the bank I in the partition J
			 */

			/** We launch the iteration of the sequences iterator with the created functors. */
			for (size_t i=0; i<itBanks.size(); i++)
			{
				size_t groupSize   = 1000;
				bool deleteSynchro = true;

				/** We fill the partitions. Each thread will read synchronously and will call FillPartitions
				 * in a synchronous way (in order to have global BanksStats correctly computed). */

			getDispatcher()->iterate(
				itBanks[i],
				FillPartitions<span, false>(
					model, _config._nb_passes, pass, _config._nb_partitions,
					_config._nb_cached_items_per_core_per_part, _progress, _bankStats,
					_tmpPartitions, *_repartitor, pInfo, _superKstorage),
				groupSize, deleteSynchro);

				/** We flush the partitions in order to be sure to have the exact number of items per partition. */
					_tmpPartitions->flush();

					/** We get a snapshot of items number in each partition. */
					vector<size_t> nbItems;
					for (size_t p=0; p<_config._nb_partitions; p++)
					{
					 nbItems.push_back ((*_tmpPartitions)[p].getNbItems()); //todo for multi count
					}

					/** We add the current number of kmers in each partition for the reached ith bank. */
					_nbKmersPerPartitionPerBank.push_back (nbItems);

				//GR: close the input bank here with call to finalize
				itBanks[i]->finalize();
			}
		}

        // force close partitions and re-open them for reading
        // may prevent crash in large multi-bank counting instance on Lustre filesystems
		if(_config._solidityKind != KMER_SOLIDITY_SUM)
		{
			string tmpStorageName = getInput()->getStr(STR_URI_OUTPUT_TMP) + "/" + System::file().getTemporaryFilename("dsk_partitions");
			setPartitions        (0); // close the partitions first, otherwise new files are opened before  closing parti from previous pass
			setPartitions        ( & (*_tmpPartitionsStorage)().getPartition<Type> ("parts"));
			
		}
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
std::vector<size_t> SortingCountAlgorithm<span>::getNbCoresList (PartiInfo<5>& pInfo)
{
    std::vector<size_t> result;

    for (size_t p=0; p<_config._nb_partitions; )
    {
        u_int64_t ram_total = 0;
        size_t i=0;
        for (i=0; i< _config._nb_partitions_in_parallel && p<_config._nb_partitions
            && (ram_total ==0  || ((ram_total+(pInfo.getNbSuperKmer(p)*getSizeofPerItem()))  <= _config._max_memory*MBYTE)) ; i++, p++)
        {
            ram_total += pInfo.getNbSuperKmer(p)*getSizeofPerItem();
        }

        result.push_back (i);
    }

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
template<size_t span>
void SortingCountAlgorithm<span>::fillSolidKmers (size_t pass, PartiInfo<5>& pInfo)
{
    TIME_INFO (getTimeInfo(), "fill_solid_kmers");

    for (size_t i=0; i<_processors.size(); i++)
    {
        /** We notify the count processor about the start of the pass. */
        _processors[i]->beginPass (pass);

        fillSolidKmers_aux (_processors[i], pass, pInfo);

        /** We notify the count processor about the end of the pass. */
        _processors[i]->endPass (pass);
    }
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
void SortingCountAlgorithm<span>::fillSolidKmers_aux (ICountProcessor<span>* processor, size_t pass, PartiInfo<5>& pInfo)
{
    DEBUG (("SortingCountAlgorithm<span>::fillSolidKmers\n"));

    /** We update the message of the progress bar. */
    _progress->setMessage (Stringify::format (progressFormat2, pass+1, _config._nb_passes));


    /** We retrieve the list of cores number for dispatching N partitions in N threads.
     *  We need to know these numbers for allocating the N maps according to the maximum allowed memory.
     */
    vector<size_t> coreList = getNbCoresList(pInfo); //uses _nb_partitions_in_parallel

    /** We need a memory allocator. We give the cores number in order to compute an extra memory
     * allocation for alignment constraints. */
    MemAllocator pool (_config._nbCores);

    size_t p = 0;
    for (size_t i=0; i<coreList.size(); i++)
    {
        vector<ICommand*> cmds;

        /** We use a vector to hold all the current CountProcessor clones. */
        vector<CountProcessor*> clones;

        size_t currentNbCores = coreList[i];
        assert (currentNbCores > 0);

        /** We correct the number of memory per map according to the max allowed memory.
         * Note that _max_memory has initially been divided by the user provided cores number. */
        u_int64_t mem = (_config._max_memory*MBYTE)/currentNbCores;

        /** We need to cache the solid kmers partitions.
         *  NOTE : it is important to save solid kmers by big chunks (ie cache size) in each partition.
         *  Indeed, if we directly iterate the solid kmers through a Partition::iterator() object,
         *  one partition is iterated after another one, which doesn't reflect the way they are in filesystem,
         *  (ie by chunks of solid kmers) which may lead to many moves into the global  h d f 5 file.
         *  One solution is to make sure that the written chunks of solid kmers are big enough: here
         *  we accept to provide at most 2% of the max memory, or chunks of 200.000 items.
         */
        size_t cacheSize = std::min ((u_int64_t)(200*1000), mem/(50*sizeof(Count)));

        DEBUG (("SortingCountAlgorithm::fillSolidKmers:  mem=%d  computing %zu partitions simultaneously , parti : ",
            mem/MBYTE, currentNbCores
        ));

        /** We build a list of 'currentNbCores' commands to be dispatched each one in one thread. */
        for (size_t j=0; j<currentNbCores; j++, p++)
        {
            ISynchronizer* synchro = System::thread().newSynchronizer();
            LOCAL (synchro);

            /** We clone the prototype count processor instance for the current 'p' kmers partition. */
            CountProcessor* processorClone = processor->clone ();

            /** We use and put the clone into a vector. */
            processorClone->use();
            clones.push_back (processorClone);

            DEBUG ((" %zu ", p));

            /* Get the memory taken by this partition if loaded for sorting */
            uint64_t memoryPartition = (pInfo.getNbSuperKmer(p)*getSizeofPerItem()); //in bytes
            DEBUG (("  (%llu  MB) ",memoryPartition/MBYTE));

            /** If we have several input banks, we may have to compute kmer solidity for each bank, which
             * can be currently done only with sorted vector. */
			bool forceVector  =   _nbKmersPerPartitionPerBank.size() > 1 && ( _config._solidityKind != KMER_SOLIDITY_SUM);


            ICommand* cmd = 0;

            //still use hash if by vector would be too large even with single part at a time
			//I thought it was not possible to have memoryPartition > _max_memory  && currentNbCores>1 , but inf fact it is possible when
			// some partitions are of size 0 (see getNbCoresList)
			if ( ((memoryPartition > mem && currentNbCores==1) || ( memoryPartition > (_config._max_memory*MBYTE) ) )  && !forceVector)
            {
                if (pool.getCapacity() != 0)  {  pool.reserve(0);  }


					cmd = new PartitionsByHashCommand<span>   (
															   processorClone, cacheSize, _progress, _fillTimeInfo,
															   pInfo, pass, p, _config._nbCores_per_partition, _config._kmerSize, pool, mem,_superKstorage
															   );
            }
            else
            {
                u_int64_t memoryPoolSize = _config._max_memory*MBYTE;

                /** In case of forcing sorted vector (multiple banks counting for instance), we may have a
                 * partition bigger than the max memory. */
                if (forceVector  &&  memoryPartition >= memoryPoolSize)
                {
                    static const int EXCEED_FACTOR = 2;

                    if (memoryPartition  < EXCEED_FACTOR*memoryPoolSize)
                    {						
                        /** We accept in this case to exceed the allowed memory. */
                        memoryPoolSize = memoryPartition;
                    }
                    else
                    {
                        bool strict = false;

                        if (strict)
                        {
                            /** We launch an exception. */
                            throw Exception ("memory issue: %lld bytes required and %lld bytes available",
                                memoryPartition, memoryPoolSize
                            );
                        }
                        else
                        {
                            unsigned long system_mem = System::info().getMemoryPhysicalTotal();
                            memoryPoolSize = memoryPartition; 

                            if (memoryPoolSize > system_mem*0.95)
                            {
                                throw Exception ("memory issue: %lld bytes required, %lld bytes set by command-line limit, %lld bytes in system memory",
                                    memoryPartition, memoryPoolSize, system_mem
                                );
                            }
                            else
                                cout << "Warning: memory was initially restricted to " << _config._max_memory << " MB, but we actually need to allocate " << memoryPoolSize / MBYTE << " MB due to a partition with " << pInfo.getNbSuperKmer(p) << " superkmers." << endl;
                        }
                    }
                }

               //if capa pool ==0, reserve max memo , pass pool to partibyvec, will be used  for vec kmers
                if (pool.getCapacity() == 0)  {  pool.reserve (memoryPoolSize); }
				else if (memoryPoolSize > pool.getCapacity()) { pool.reserve(0); pool.reserve (memoryPoolSize); }

                /** Recall that we got the following matrix in _nbKmersPerPartitionPerBank
                 *
                 *           part0  part1  part2 ... partJ
                 *   bank0    xxx    xxx    xxx       xxx
                 *   bank1    xxx    xxx    xxx       xxx
                 *    ...
                 *   bankI    xxx    xxx    xxx       xxx
                 *
                 *   Now, for the current partition p, we want the number of items found for each bank.
                 *
                 *              bank0   bank1   ...   bankI
                 *   offsets :   xxx     xxx           xxx
                 */
                vector<size_t> nbItemsPerBankPerPart;
                if ( _config._solidityKind != KMER_SOLIDITY_SUM)
                {
                    for (size_t i=0; i<_nbKmersPerPartitionPerBank.size(); i++)
                    {
                        nbItemsPerBankPerPart.push_back (_nbKmersPerPartitionPerBank[i][p] - (i==0 ? 0 : _nbKmersPerPartitionPerBank[i-1][p]) );
                    }
                }

				if ( _config._solidityKind == KMER_SOLIDITY_SUM)
				{
					cmd = new PartitionsByVectorCommand<span> (
															   processorClone, cacheSize, _progress, _fillTimeInfo,
															   pInfo, pass, p, _config._nbCores_per_partition, _config._kmerSize, pool, nbItemsPerBankPerPart,_superKstorage
															   );
				}
				else
				{
					cmd = new PartitionsByVectorCommand_multibank<span> (
															   (*_tmpPartitions)[p], processorClone, cacheSize, _progress, _fillTimeInfo,
															   pInfo, pass, p, _config._nbCores_per_partition, _config._kmerSize, pool, nbItemsPerBankPerPart
															   );
				}

            }

            cmds.push_back (cmd);

        } /* end of for (size_t j=0; j<currentNbCores... */

        DEBUG (("\n"));

        /** We launch the commands through a dispatcher. */
        getDispatcher()->dispatchCommands (cmds, 0);

        /** The N CountProcessor clones should have done their job during the 'dispatchCommands'
         * We can send a notification about it and get rid of them. */
        processor->finishClones (clones);
        for (size_t i=0; i<clones.size(); i++)  { clones[i]->forget(); }  clones.clear();

        // free internal memory of pool here
        pool.free_all();
    }
	
	
	if(_config._solidityKind == KMER_SOLIDITY_SUM)
		_superKstorage->closeFiles();

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
Partition<typename SortingCountAlgorithm<span>::Count>* SortingCountAlgorithm<span>::getSolidCounts  ()
{
    /** We look in the count processor a potential CountProcessorDump instance. */
    for (size_t i=0; i<_processors.size(); i++)
    {
        CountProcessorDump<span>* p = _processors[i]->template get <CountProcessorDump<span> > ();
        if (p != 0)  {  return p->getSolidCounts();  }
    }

    throw Exception ("SortingCountAlgorithm not configured with a CountProcessorDump instance");
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
struct Count2TypeAdaptor  {  typename Kmer<span>::Type& operator() (typename Kmer<span>::Count& c)  { return c.value; }  };

template<size_t span>
Iterable<typename SortingCountAlgorithm<span>::Type>* SortingCountAlgorithm<span>::getSolidKmers ()
{
    Iterable<Type>* result = 0;

    tools::storage::impl::Partition<Count>* counts = getSolidCounts();

    if (counts != 0) { result = new IterableAdaptor<Count,Type,Count2TypeAdaptor<span> > (*counts); }

    return result;
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
