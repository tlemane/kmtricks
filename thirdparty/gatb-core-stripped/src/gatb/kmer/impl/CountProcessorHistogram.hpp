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

#ifndef _COUNT_PROCESSOR_HISTOGRAM_HPP_
#define _COUNT_PROCESSOR_HISTOGRAM_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/CountProcessorAbstract.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/misc/impl/Histogram.hpp>
#include <cstdarg>
#include <gatb/system/impl/System.hpp>

using namespace gatb::core::system;
using namespace gatb::core::system::impl;


/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

/** The CountProcessorHistogram implementation collects information about the
 * kmers distribution. It feeds a IHistogram instance during the 'process' method.
 *
 * At the end of the algorithm, it provides information such as the best cutoff
 * (ie a "good" abundance min parameter).
 */
template<size_t span=KMER_DEFAULT_SPAN>
class CountProcessorHistogram : public CountProcessorAbstract<span>
{
public:

    /** Shortcuts. */
    typedef typename Kmer<span>::Type Type;

    /** Constructor. */
    CountProcessorHistogram (
        tools::storage::impl::Group* group = 0,
        size_t histoMax                    = 10000,
        size_t min_auto_threshold          = 3,
		bool   histo2Dmode                 = false,
		bool   histo1Dmode                 = false,
		std::string histo2Dfilename = "histo2Dfile",
		std::string histo1Dfilename = "histo1Dfile"
    )
        : _group(group), _histogram(0), _min_auto_threshold(min_auto_threshold),_histo2Dmode(histo2Dmode),_histo2Dfilename(histo2Dfilename),_histo1Dmode(histo1Dmode),_histo1Dfilename(histo1Dfilename)
    {

        setHistogram (new tools::misc::impl::Histogram (histoMax));
		_synchro = System::thread().newSynchronizer(); // synchro that will be passed down to histogramCache in each clone

    }

    /** Constructor. */
    CountProcessorHistogram (
        tools::storage::impl::Group* group,
        tools::misc::IHistogram* histogram,
        size_t min_auto_threshold = 3,
		bool   histo2Dmode = false,
		bool   histo1Dmode = false,
		std::string histo2Dfilename = "histo2Dfile",
	    std::string histo1Dfilename = "histo1Dfile"
    )
        : _group(group), _histogram(0), _min_auto_threshold(min_auto_threshold),_histo2Dmode(histo2Dmode),_histo2Dfilename(histo2Dfilename),_histo1Dmode(histo1Dmode),_histo1Dfilename(histo1Dfilename)
    {

        setHistogram (histogram);
		_synchro = 0 ; // no need for another synchro object in the clones
    }

    /** Destructor. */
    virtual ~CountProcessorHistogram ()
    {
		delete _synchro;
        setHistogram(0);
    }

    /********************************************************************/
    /*   METHODS CALLED ON THE PROTOTYPE INSTANCE (in the main thread). */
    /********************************************************************/

    /** \copydoc ICountProcessor<span>::end */
    void end ()
    {
        using namespace tools::math;

        /** compute auto cutoff **/
        _histogram->compute_threshold (_min_auto_threshold);

		if(_histo2Dmode)
		{
			FILE * histo2Dfile = fopen (_histo2Dfilename.c_str(),"w");
			//output 2D histogram now
			//printf("output 2D histo gram to file %s \n",_histo2Dfilename.c_str());
			
			for(int ii=0; ii<= _histogram->getLength(); ii++)
			{
				fprintf(histo2Dfile,"%5i:\t",ii);
				for(int jj=0; jj<= _histogram->getLength2(); jj++)
				{
					fprintf(histo2Dfile,"\t%6lli", _histogram->get2D(ii,jj));
				}
				fprintf(histo2Dfile,"\n");
			}
			
			fclose(histo2Dfile);
		}
		
		if(_histo1Dmode)
		{
			FILE * histo1Dfile = fopen (_histo1Dfilename.c_str(),"w");
			
			//output 1D histogram now
			for(int ii=1; ii<= _histogram->getLength(); ii++)
			{
				fprintf(histo1Dfile,"%i\t%lli",ii,_histogram->get(ii));
				fprintf(histo1Dfile,"\n");
			}
			
			fclose(histo1Dfile);
		}
		
		
        if (_group != 0)
        {
            /** We save the histogram if any. */
            _histogram->save (*_group);

            /** store auto cutoff and corresponding number of solid kmers **/
            tools::collections::Collection<NativeInt64>& storecutoff = _group->getCollection<NativeInt64>("cutoff") ;
            storecutoff.insert(_histogram->get_solid_cutoff());
            storecutoff.flush();

            tools::collections::Collection<NativeInt64>& storesolids = _group->getCollection<NativeInt64>("nbsolidsforcutoff") ;
            storesolids.insert(_histogram->get_nbsolids_auto());
            storesolids.flush();
        }
    }

    /** \copydoc ICountProcessor<span>::clone */
    CountProcessorAbstract<span>* clone ()
    {
        /** We encapsulate the histogram with a cache. */
        return new CountProcessorHistogram (_group, new gatb::core::tools::misc::impl::HistogramCache (_histogram,_synchro),  _min_auto_threshold, _histo2Dmode,_histo1Dmode, _histo2Dfilename,_histo1Dfilename);
    }

    /********************************************************************/
    /*   METHODS CALLED ON ONE CLONED INSTANCE (in a separate thread).  */
    /********************************************************************/

    /** \copydoc ICountProcessor<span>::process */
    bool process (size_t partId, const Type& kmer, const CountVector& count, CountNumber sum)
    {
        _histogram->inc (sum);
		
		if(_histo2Dmode)
		{
			CountNumber sumreads = sum - count[0];
			_histogram->inc2D (sumreads,count[0]);
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

        result.add (0, "histogram");
        result.add (1, "cutoff",            "%ld",  _histogram->get_solid_cutoff());
        result.add (1, "nb_ge_cutoff",      "%ld",  _histogram->get_nbsolids_auto());
		result.add (1, "ratio_weak_volume",      "%.2f",  _histogram->get_ratio_weak());
		
		
        // result->add (1, "percent_ge_cutoff", "%.1f", nbSolids > 0 ? 100.0 * (double)_histogram->get_nbsolids_auto() / (double)_bankStats.kmersNbValid : 0);
        result.add (1, "first_peak",         "%ld",  _histogram->get_first_peak());

        // double N = ((double)_histogram->get_first_peak() * _bankStats.getSeqMean()) / (_bankStats.getSeqMean() - _kmerSize + 1);
        // if (N > 0)  {  getInfo()->add (3, "genome_size_estimate", "%.0f",  (double)_bankStats.sequencesTotalLength / N);  }

        return result;
    }

    /** Get the histogram.
     * \return the histogram instance. */
    gatb::core::tools::misc::IHistogram* getHistogram() { return _histogram; }

private:


	system::ISynchronizer* _synchro;

    tools::storage::impl::Group* _group;

    gatb::core::tools::misc::IHistogram* _histogram;
    void setHistogram (gatb::core::tools::misc::IHistogram* histogram)  { SP_SETATTR(histogram); }

    size_t _min_auto_threshold;
	bool _histo2Dmode;
	std::string _histo2Dfilename;
	bool _histo1Dmode;
	std::string _histo1Dfilename;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _COUNT_PROCESSOR_HISTOGRAM_HPP_ */
