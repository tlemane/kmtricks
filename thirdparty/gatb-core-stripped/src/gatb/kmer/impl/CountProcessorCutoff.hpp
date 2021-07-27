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

#ifndef _COUNT_PROCESSOR_CUTOFF_HPP_
#define _COUNT_PROCESSOR_CUTOFF_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/CountProcessorAbstract.hpp>
#include <gatb/kmer/impl/CountProcessorHistogram.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/misc/impl/Histogram.hpp>
#include <cstdarg>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

/** The CountProcessorCutoff implementation collects information about the
 * kmers distribution. It feeds a IHistogram instance during the 'process' method.
 *
 * At the end of the algorithm, it provides information such as the best cutoff
 * (ie a "good" abundance min parameter).
 */
template<size_t span=KMER_DEFAULT_SPAN>
class CountProcessorCutoff : public CountProcessorAbstract<span>
{
public:

    /** Shortcuts. */
    typedef typename Kmer<span>::Type Type;

    /** Constructor for the prototype instance. */
    CountProcessorCutoff (size_t nbBanks) : CountProcessorAbstract<span>("cutoff")
    {
        /** We create one histogram instance per bank. */
        for (size_t i=0; i<nbBanks; i++)  {  _histogramProcessors.push_back (new CountProcessorHistogram<span> ());  }
    }

    /** Constructor for the cloned instances. */
    CountProcessorCutoff (const std::vector <CountProcessorHistogram<span>* >& histogramProcessors)
        : CountProcessorAbstract<span>("cutoff"), _histogramProcessors(histogramProcessors) {}

    /** Destructor. */
    virtual ~CountProcessorCutoff ()
    {
        for (size_t i=0; i<_histogramProcessors.size(); i++)  {  delete _histogramProcessors[i];  }
    }

    /********************************************************************/
    /*   METHODS CALLED ON THE PROTOTYPE INSTANCE (in the main thread). */
    /********************************************************************/

    /** \copydoc ICountProcessor<span>::clone */
    CountProcessorAbstract<span>* clone ()
    {
        // We clone each CountProcessorHistogram
        std::vector <CountProcessorHistogram<span>* > clones;
        for (size_t i=0; i<_histogramProcessors.size(); i++)  {  clones.push_back ((CountProcessorHistogram<span>*)_histogramProcessors[i]->clone());  }

        // We return an instance of our custom processor with the CountProcessorHistogram clones
        return new CountProcessorCutoff (clones);
    }

    /** \copydoc ICountProcessor<span>::end */
    void endPass (size_t passId)
    {
        size_t _min_auto_threshold = 3;

        _cutoffs.clear();

        for (size_t i=0; i<_histogramProcessors.size(); i++)
        {
            /** Shortcut. */
            gatb::core::tools::misc::IHistogram* histogram = _histogramProcessors[i]->getHistogram();

            histogram->compute_threshold (_min_auto_threshold);

            _cutoffs.push_back (histogram->get_solid_cutoff());
        }
    }

    /********************************************************************/
    /*   METHODS CALLED ON ONE CLONED INSTANCE (in a separate thread).  */
    /********************************************************************/

    /** \copydoc ICountProcessor<span>::process */
    bool process (size_t partId, const Type& kmer, const CountVector& count, CountNumber sum)
    {
        // In case we have only one bank, we use the sum by convention. Note that it can still
        // be the sum of several bank counts (here the _histogramProcessors.size() value may not reflect
        // the actual number of banks processed, see SortingCountAlgorithm<span>::getDefaultProcessorVector)

        if (_histogramProcessors.size()==1)
        {
            sum=0; for (size_t i=0; i<count.size(); i++)  { sum += count[i]; }
            _histogramProcessors[0]->process (partId, kmer, count, sum);
        }
        else
        {
            for (size_t i=0; i<_histogramProcessors.size(); i++)  {  _histogramProcessors[i]->process (partId, kmer, count, count[i]);  }
        }
        return true;
    }

    /*****************************************************************/
    /*                          MISCELLANEOUS.                       */
    /*****************************************************************/

    /** \copydoc ICountProcessor<span>::getProperties */
    tools::misc::impl::Properties getProperties() const
    {
        tools::misc::impl::Properties result;

        std::stringstream ss;  for (size_t i=0; i<_cutoffs.size(); i++)  { ss << _cutoffs[i] << " "; }
        result.add (0, "values", "%s", ss.str().c_str());

        return result;
    }

    /** */
    std::vector<CountNumber> getCutoffs()  {  return _cutoffs;  }

private:

    std::vector<CountProcessorHistogram<span>* > _histogramProcessors;

    std::vector<CountNumber> _cutoffs;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _COUNT_PROCESSOR_CUTOFF_HPP_ */
