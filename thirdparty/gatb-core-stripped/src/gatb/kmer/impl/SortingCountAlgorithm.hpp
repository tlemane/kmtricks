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

/** \file SortingCountAlgorithm.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Counting kmers from a set of sequences
 */

#ifndef _SORTING_COUNT_ALGORITHM_HPP_
#define _SORTING_COUNT_ALGORITHM_HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/bank/api/IBank.hpp>
#include <gatb/kmer/api/ICountProcessor.hpp>
#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/BankKmers.hpp>
#include <gatb/kmer/impl/Configuration.hpp>
#include <gatb/kmer/impl/PartiInfo.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>
#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for genomic databases management. */
namespace kmer      {
/** \brief Implementation for genomic databases management. */
namespace impl      {
/********************************************************************************/

/** \brief Class performing the kmer counting (also known as 'DSK')
 *
 * This class does the real job of counting the kmers from a reads database.
 *
 * This is a template class whose template argument is the kind of integer used for
 * kmers (integers on 64 bits, 128 bits, etc...)
 *
 * We define some template instantiations of this SortingCountAlgorithm; such an instantiation
 * does the real job of kmers counting. By defining several instantiations, we allow
 * to choose dynamically the correct class according to the user choice for kmer size
 * (remember that initial Minia version had to be re-compiled for different kmer size).
 *
 * Actually, this class is mainly used in the debruijn::impl::Graph class as a first step for
 * the de Bruijn graph creation.
 */
template<size_t span=KMER_DEFAULT_SPAN>
class SortingCountAlgorithm : public gatb::core::tools::misc::impl::Algorithm
{
public:

	/* Shortcuts. */
    typedef typename Kmer<span>::Type                                       Type;
    typedef typename Kmer<span>::ModelDirect                            ModelDirect;
    typedef typename Kmer<span>::ModelCanonical                         ModelCanonical;
#ifdef NONCANONICAL
    typedef typename Kmer<span>::template ModelMinimizer <ModelDirect>   Model;
#else
    typedef typename Kmer<span>::template ModelMinimizer <ModelCanonical>   Model;
#endif
    typedef typename Kmer<span>::Count                                      Count;
    typedef ICountProcessor<span> CountProcessor;

    /** Constructor. Can be used as default constructor in no parameters are provided
     * \param[in] params : parameters to be used for configuring the algorithm
     */
    SortingCountAlgorithm (tools::misc::IProperties* params = 0);

    /** Constructor.
     * \param[in] bank : input bank from which solid kmers are counted
     * \param[in] params : parameters to be used for configuring the algorithm
     */
    SortingCountAlgorithm (
        gatb::core::bank::IBank*  bank,
        tools::misc::IProperties* params
    );

    /** Constructor.
     * \param[in] bank : input bank from which solid kmers are counted
     * \param[in] config : configuration information
     * \param[in] repartitor : hash function for minimizers
     * \param[in] processor : object that processes counts
     */
    SortingCountAlgorithm (
        gatb::core::bank::IBank*     bank,
        const Configuration&         config,
        Repartitor*                  repartitor,
        std::vector<CountProcessor*> processors,
		tools::misc::IProperties* params

    );

    /** Destructor */
    virtual ~SortingCountAlgorithm ();

    /** operator=
     * \param[in] s : object to be copied. */
    SortingCountAlgorithm& operator= (const SortingCountAlgorithm& s);

    /** Get an option parser for kmers counting parameters. Dynamic allocation, so must be released when no more used.
     * \param[in] mandatory : tells whether an argument has to be mandatory
     * \return an instance of IOptionsParser. */
    static tools::misc::IOptionsParser* getOptionsParser (bool mandatory=true);

    /** Get the default values defined in the default option parser.
     * \return default properties. */
    static tools::misc::IProperties* getDefaultProperties ();

    /** Creates a default CountProcessor instance (ie. the default one used by DSK)
     * \param[in] params : used for configuring the processor
     * \param[in] dskStorage : storage for dumping [kmer,count] couples
     * \param[in] otherStorage : used for histogram for instance
     * \return a CountProcessor instance
     */
    static CountProcessor* getDefaultProcessor (
        tools::misc::IProperties*       params,
        tools::storage::impl::Storage*  dskStorage,
        tools::storage::impl::Storage*  otherStorage = 0
    );

    /** Creates a vector holding the default CountProcessor configuration
     * \param[in] params : used for configuring the processor
     * \param[in] dskStorage : storage for dumping [kmer,count] couples
     * \param[in] otherStorage : used for histogram for instance
     * \return a vector of CountProcessor instances
     */
    static std::vector<ICountProcessor<span>*> getDefaultProcessorVector (
        Configuration&                  config,
        tools::misc::IProperties*       params,
        tools::storage::impl::Storage*  dskStorage,
        tools::storage::impl::Storage*  otherStorage = 0
    );

    /** Process the kmers counting. It is mainly composed of a loop over the passes, and for each pass :
     *      1) we build the partition files then
     *      2) we fill the solid kmers file from the partitions.
     */
    void  execute ();

    /** Get the number of count processors associated to the object.
     * \return number of processors. */
    size_t getProcessorNumber() const { return _processors.size(); }

    /** Getter for the CountProcessor to be used by the algorithm.
     * \param[in] idx : index of the processor to be retrieved
     * \return the processor
     */
    CountProcessor* getProcessor (size_t idx)  { return _processors[idx]; }

    /** Setter for the CountProcessor to be used by the algorithm.
     * \param[in] processor : the count processor to be used.
     */
    void addProcessor (CountProcessor* processor)  { processor->use(); _processors.push_back (processor); }

    /** Get the iterable over the computed solid kmers.
     * \return the solid kmers iterable. */
    tools::storage::impl::Partition<Count>* getSolidCounts ();

    /** Get the iterable over the computed solid kmers.
     * \return the solid kmers iterable. */
    tools::collections::Iterable<Type>* getSolidKmers   ();

    /** Get the configuration object for the algorithm.
     * \return the Configuration object. */
    const kmer::impl::Configuration& getConfig() const { return _config; }

    /* Get the storage instance (if any).
     * \return the Storage instance. */
    tools::storage::impl::Storage* getStorage () { return _storage; }

    /** Get the repartitor instance, ie. the hash function built on minimizer information.
     * \return the Repartitor instance
     */
    Repartitor* getRepartitor() { return _repartitor; }

private:

    /** Configuration of the objects used by the algorithm. */
    void configure ();

    /** Fill partition files (for a given pass) from a sequence iterator.
     * \param[in] pass  : current pass whose value is used for choosing the partition file
     * \param[in] itSeq : sequences iterator whose sequence are cut into kmers to be split.
     */
    void fillPartitions (size_t pass, gatb::core::tools::dp::Iterator<gatb::core::bank::Sequence>* itSeq, PartiInfo<5>& pInfo);

    /** Fill the solid kmers bag from the partition files (one partition after another one).
     * \param[in] solidKmers : bag to put the solid kmers into.
     */
    void fillSolidKmers (size_t pass, PartiInfo<5>& pInfo);

    /** Fill the solid kmers bag from the partition files (one partition after another one).
     * \param[in] solidKmers : bag to put the solid kmers into.
     */
    void fillSolidKmers_aux (ICountProcessor<span>* processor, size_t pass, PartiInfo<5>& pInfo);

    /** */
    std::vector <size_t> getNbCoresList (PartiInfo<5>& pInfo);

    /** Handle on the configuration information. */
    kmer::impl::Configuration _config;

    /** Handle on the input bank. */
    gatb::core::bank::IBank* _bank;
    void setBank (gatb::core::bank::IBank* bank)  { SP_SETATTR(bank); }

    /** Handle on the mininimizers hash function. */
    Repartitor* _repartitor;
    void setRepartitor (Repartitor* repartitor)  { SP_SETATTR(repartitor); }

    /** Handle on the count processor object. */
    std::vector<CountProcessor*> _processors;

	
    /** Handle on the progress information. */
    gatb::core::tools::dp::IteratorListener* _progress;
    void setProgress (gatb::core::tools::dp::IteratorListener* progress)  { SP_SETATTR(progress); }

    /** Temporary partitions management. */
    tools::storage::impl::Storage* _tmpPartitionsStorage;
    void setPartitionsStorage (tools::storage::impl::Storage* tmpPartitionsStorage)  {  SP_SETATTR(tmpPartitionsStorage);  }

    /** Temporary partitions management. */
    tools::storage::impl::Partition<Type>* _tmpPartitions;
    void setPartitions (tools::storage::impl::Partition<Type>* tmpPartitions)  {  SP_SETATTR(tmpPartitions);  }

    /** Get the memory size (in bytes) to be used by each item.
     * IMPORTANT : we may have to count both the size of Type and the size for the bank id. */
    int getSizeofPerItem () const { return Type::getSize()/8 + ((_nbKmersPerPartitionPerBank.size()>1 && _config._solidityKind != tools::misc::KMER_SOLIDITY_SUM) ? sizeof(bank::BankIdType) : 0); }

    tools::misc::impl::TimeInfo _fillTimeInfo;

    BankStats _bankStats;

    std::vector <std::vector<size_t> > _nbKmersPerPartitionPerBank;

    tools::storage::impl::StorageMode_e _storage_type;
    tools::storage::impl::Storage* _storage;
    void setStorage (tools::storage::impl::Storage* storage)  { SP_SETATTR(storage); }
	
	
	//superkmer efficient storage
	tools::storage::impl::SuperKmerBinFiles* _superKstorage;
	std::string _tmpStorageName_superK;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _SORTING_COUNT_ALGORITHM_HPP_ */

