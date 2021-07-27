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

#ifndef _GATB_CORE_KMER_ICOUNT_PROCESSOR_HPP_
#define _GATB_CORE_KMER_ICOUNT_PROCESSOR_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/Configuration.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
/********************************************************************************/

/** \brief Interface that uses kmer counting information
 *
 * This interface is mainly an Observer that listens to data produced by the sorting
 * count algorithm. Such an information is made of a kmer an the number of occurrences
 * of this kmers in each bank provided to the algorithm.
 *
 * Through this interface, it becomes easy to plug specific listeners that can do
 * different things on the [kmer,counts] information. There is a default implementation
 * of the ICountProcessor interface that does the historical job of DSK:
 *      1) building an histogram
 *      2) filtering out kmers with too low coverage
 *      3) saving on disk kmers having big enough coverage
 *
 * Such an instance can be associated to the SortingCountAlgorithm instance with the
 * SortingCountAlgorithm::setProcessor method; this instance will be called 'prototype
 * instance'.
 *
 * From an execution point of view, one instance of ICountProcessor is created (with method
 * 'clone') for counting the kmers of one specific partition. If N cores are used, it means
 * that N instances of ICountProcessor will be cloned from the so called 'prototype'
 * instance (ie. the instance associated to the SortingCountAlgorithm instance). Each
 * clone processes its partition in one specific thread.
 *
 * While processing a partition, a cloned ICountProcessor instance is called
 * via its 'process' method: this is here that the information [kmer,counts] is
 * provided to the ICountProcessor clone, and accordingly to the actual implementation
 * class of the ICountProcessor interface, different processings can be done.
 *
 * When all the clones have finished their job (in their own thread), the prototype
 * instance is called (in the main thread) via the 'finishClones' method, where
 * the prototype instance has access to the N clones before they are deleted. It allows
 * for instance to gather in the prototype instance the information collected by the clones
 * during their processing.
 *
 * From a global point of view, the interface is made of three parts :
 *      1) methods called on the prototype instance in the context of the main thread
 *      2) methods called on a cloned instance in the context of specific threads
 *      3) all other methods
 *
 * The following figure shows how ICountProcessor interacts with other classes, and in particular
 * with SortingCountAlgorithm. One can also see the multithreading context, with the main thread
 * creating clones and with clones processing their job in specific threads.
 *
 * \image html ICountProcessor.png "Usage and life cycle of ICountProcessor in the context of SortingCountAlgorithm"
 *
 * Examples of ICountProcessor implementors :
 *      1) CountProcessorHistogram   : collect kmers distribution information
 *      2) CountProcessorSolidity... : check whether a kmer is solid or not
 *      3) CountProcessorDump        : dump kmer count information in file system
 *      4) CountProcessorChain       : list of linked ICountProcessor instances
 *
 * The CountProcessorChain implementation allows to link several instances of
 * ICountProcessor. When such an instance is called via 'process', the first item
 * of the list is called via 'process'; if it returns true, the next item in the list
 * is called and so on; if it returns false, the chain is stopped. This class is used
 * for the definition of the "DSK" count processor (histogram -> solidity -> dump)
 */
template<size_t span>
class ICountProcessor : public system::SmartPointer
{
public:

    /** Shortcuts. */
    typedef typename kmer::impl::Kmer<span>::Type Type;

    /********************************************************************/
    /*   METHODS CALLED ON THE PROTOTYPE INSTANCE (in the main thread). */
    /********************************************************************/

    /** Called just before the mainloop of SortingCountAlgorithm.
     * \param[in] config : configuration of the SortingCountAlgorithm. */
    virtual void begin (const kmer::impl::Configuration& config) = 0;

    /** Called just after the mainloop of SortingCountAlgorithm. */
    virtual void end   () = 0;

    /** Called just before starting a pass.
     * \param[in] passId: index of the pass to begin */
    virtual void beginPass (size_t passId) = 0;

    /** Called just after the end of a pass. */
    virtual void endPass   (size_t passId) = 0;

    /** Clone the instance.
     * An instance can be cloned N times in order to use the cloned instance in one thread.
     * \return the cloned instance. */
    virtual ICountProcessor* clone () = 0;

    /** Called when N partitions have been processed through N clones. This should be the last
     * time these clones are available before being deleted. It can be the opportunity to the
     * prototype instance to gather information from the clones.
     * \param[in] clones : the N cloned instances
     */
    virtual void finishClones (std::vector<ICountProcessor<span>*>& clones) = 0;

    /********************************************************************/
    /*   METHODS CALLED ON ONE CLONED INSTANCE (in a separate thread).  */
    /********************************************************************/

    /** Called at the beginning of a new kmers partition processing.
     * \param[in] passId : index of the current pass in the SortingCountAlgorithm.
     * \param[in] passId : index of the current kmers partition in the SortingCountAlgorithm.
     * \param[in] cacheSize : memory size used for the current kmers partition
     * \param[in] name : class name of the child PartitionsCommand class.
     */
    virtual void beginPart (size_t passId, size_t partId, size_t cacheSize, const char* name) = 0;

    /** Called at the end of a new kmers partition processing.
     * \param[in] passId : index of the current pass in the SortingCountAlgorithm.
     * \param[in] passId : index of the current kmers partition in the SortingCountAlgorithm.
     */
    virtual void endPart   (size_t passId, size_t partId) = 0;

    /** Notification that a [kmer,counts] is available and can be handled by the count processor.
     * \param[in] partId : index of the current partition
     * \param[in] kmer : kmer for which we are receiving counts
     * \param[in] count : vector of counts of the kmer, one count per bank
     * \param[in] sum : sum of the occurrences for all bank.
     */
    virtual bool process (size_t partId, const Type& kmer, const CountVector& count, CountNumber sum=0) = 0;

    /*****************************************************************/
    /*                          MISCELLANEOUS.                       */
    /*****************************************************************/

    /** Get a name for the count processor.
     * \return the count processor name. */
    virtual std::string getName() const = 0;

    /** Set a name for the count processor.
     * \param[in] name : the count processor name. */
    virtual void setName (const std::string& name) = 0;

    /** Get some properties about the count processor.
     * \return properties. */
    virtual tools::misc::impl::Properties getProperties() const = 0;

    /** Get a vector of instances in case of the current object is a composite.
     * \return a vector of ICountProcessor instance. */
    virtual std::vector<ICountProcessor*> getInstances () const = 0;

    /** Try to get an instance of a specific type within the current object.
     * \return a T pointer to the instance if found, 0 otherwise. */
    template<typename T> T* get () const
    {
        std::vector<ICountProcessor*> v = this->getInstances();
        for (size_t i=0; i<v.size(); i++)  {  if (T* object = dynamic_cast<T*> (v[i]))  { return object; }  }
        return (T*)0;
    }
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_ICOUNT_PROCESSOR_HPP_ */
