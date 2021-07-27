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

#ifndef _COUNT_PROCESSOR_DUMP_HPP_
#define _COUNT_PROCESSOR_DUMP_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/CountProcessorAbstract.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

/** The CountProcessorDump implementation dumps kmers on file system through
 * a Partition object.
 *
 * The number of partitions is not known at creation and is received through a
 * call to the 'begin' method; there, the number of partitions is used to
 * configure the Partition object.
 *
 * A clone instance is dedicated to one partition, which is configured when the
 * method 'beginPart' is called on a cloned instance.
 *
 * The CountProcessorDump implementation is likely to be used in a CountProcessorChain,
 * like this : solidity -> dump.  It allows to dump on file system only solid kmers.
 */
template<size_t span=KMER_DEFAULT_SPAN>
class CountProcessorDump : public CountProcessorAbstract<span>
{
public:

    /** Shortcuts. */
    typedef typename Kmer<span>::Count Count;
    typedef typename Kmer<span>::Type Type;

    /** Constructor */
    CountProcessorDump (
        tools::storage::impl::Group&            group,
        size_t                                  kmerSize,
        system::ISynchronizer*                  synchronizer = 0,
        tools::storage::impl::Partition<Count>* solidCounts  = 0,
        size_t                                  nbPartsPerPass = 0
    )
        : _group(group), _kmerSize(kmerSize), _nbPartsPerPass(nbPartsPerPass), _synchronizer(0), _solidCounts(0), _solidKmers(0)
    {
        setSolidCounts  (solidCounts);
        setSynchronizer (synchronizer);
    }

    /** Destructor */
    virtual ~CountProcessorDump ()
    {
        setSolidCounts(0);
        setSolidKmers (0);
    }

    /********************************************************************/
    /*   METHODS CALLED ON THE PROTOTYPE INSTANCE (in the main thread). */
    /********************************************************************/

    /** \copydoc ICountProcessor<span>::begin */
    void begin (const Configuration& config)
    {
        /** We remember the number of partitions for one pass. */
        _nbPartsPerPass = config._nb_partitions;

        /** We compute the number of partitions. */
        size_t nbTotalPartitions = config._nb_partitions * config._nb_passes;

        // We create the partition into the dsk group
        setSolidCounts (& _group.getPartition<Count> ("solid", nbTotalPartitions));

        /** We save (as metadata) some information. */
        _group.addProperty ("kmer_size", tools::misc::impl::Stringify::format("%d", _kmerSize));
    }

    /** \copydoc ICountProcessor<span>::clones */
    CountProcessorAbstract<span>* clone ()
    {
        /** Note : we share the synchronizer for all the clones. */
        return new CountProcessorDump (_group, _kmerSize, _synchronizer, _solidCounts, _nbPartsPerPass);
    }

    /** \copydoc ICountProcessor<span>::finishClones */
    void finishClones (std::vector<ICountProcessor<span>*>& clones)
    {
        for (size_t i=0; i<clones.size(); i++)
        {
            /** We have to recover type information. */
            if (CountProcessorDump* clone = dynamic_cast<CountProcessorDump*> (clones[i]))
            {
                for (std::map<std::string,size_t>::iterator it = clone->_namesOccur.begin(); it != clone->_namesOccur.end(); ++it)
                {
                    this->_namesOccur[it->first] += it->second;
                }
            }
        }
    }

    /********************************************************************/
    /*   METHODS CALLED ON ONE CLONED INSTANCE (in a separate thread).  */
    /********************************************************************/

    /** \copydoc ICountProcessor<span>::beginPart */
    void beginPart (size_t passId, size_t partId, size_t cacheSize, const char* name)
    {
        /** We get the actual partition idx in function of the current partition AND pass identifiers. */
        size_t actualPartId = partId + (passId * _nbPartsPerPass);

        /** We get a handle on the current solid bag (with a cache). */
        setSolidKmers (new tools::collections::impl::BagCache<Count> (& (*_solidCounts)[actualPartId], cacheSize, _synchronizer));

        /** We update some stats (want to know how many "hash" or "vector" partitions we use). */
        _namesOccur[name] ++;
    }

    /** \copydoc ICountProcessor<span>::endPart */
    void endPart (size_t passId, size_t partId)
    {
        /** We flush the current collection for the partition just processed. */
        _solidKmers->flush();
    }

    /** \copydoc ICountProcessor<span>::process */
    bool process (size_t partId, const Type& kmer, const CountVector& count, CountNumber sum)
    {
        this->_solidKmers->insert (Count(kmer,sum));
        return true;
    }

    /*****************************************************************/
    /*                          MISCELLANEOUS.                       */
    /*****************************************************************/

    /** \copydoc ICountProcessor<span>::getProperties */
    tools::misc::impl::Properties getProperties() const
    {
        tools::misc::impl::Properties result;

        size_t smallestPartition = ~0;
        size_t biggestPartition  = 0;
        for (size_t i=0; i<_solidCounts->size(); i++)
        {
            size_t currentNb = (*_solidCounts)[i].getNbItems();
            smallestPartition = std::min (smallestPartition, currentNb);
            biggestPartition  = std::max (biggestPartition,  currentNb);
        }

        result.add (0, "partitions");
        result.add (1, "nb_partitions", "%ld", _solidCounts->size());
        result.add (1, "nb_items",      "%ld", _solidCounts->getNbItems());
        result.add (1, "part_biggest",  "%ld", biggestPartition);
        result.add (1, "part_smallest", "%ld", smallestPartition);

        if (_solidCounts->size())
        {
            result.add (1, "part_mean", "%.1f", (double)_solidCounts->getNbItems() / (double)_solidCounts->size());
        }

        result.add (1, "kind");
        for (std::map<std::string,size_t>::const_iterator it = _namesOccur.begin(); it != _namesOccur.end(); ++it)
        {
            result.add (2, it->first.c_str(), "%d", it->second);
        }

        return result;
    }

    /** Get the partition of counts.
     * \return the Partition<Count> instance */
    tools::storage::impl::Partition<Count>* getSolidCounts () { return _solidCounts; }

    /** Get the number of items.
     * \return the total number of items in the partition. */
    u_int64_t getNbItems ()  { return _solidCounts ? _solidCounts->getNbItems() : 0; }

private:

    tools::storage::impl::Group& _group;

    size_t _kmerSize;

    size_t _nbPartsPerPass;

    system::ISynchronizer* _synchronizer;
    void setSynchronizer (system::ISynchronizer* synchronizer)  { SP_SETATTR(synchronizer); }

    tools::storage::impl::Partition<Count>* _solidCounts;
    void setSolidCounts (tools::storage::impl::Partition<Count>* solidCounts) { SP_SETATTR(solidCounts); }

    tools::collections::Bag<Count>* _solidKmers;
    void setSolidKmers (tools::collections::Bag<Count>* solidKmers)  {  SP_SETATTR(solidKmers);  }

    std::map<std::string,size_t> _namesOccur;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _COUNT_PROCESSOR_DUMP_HPP_ */
