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

#include <gatb/tools/misc/impl/Histogram.hpp>

#include <stdarg.h>
#include <stdio.h>

#define DEBUG(a)  //printf a

using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;

/********************************************************************************/
namespace gatb {  namespace core { namespace tools {  namespace misc {  namespace impl {
/********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Histogram::save (tools::storage::impl::Group& group)
{
    DEBUG (("Histogram::save  size=%ld\n", _length+1));

    tools::collections::Collection<Entry>& collection = group.getCollection<Entry> ("histogram");

    size_t offset = 1;
    collection.insert (_histogram + offset, (_length+1) - offset);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Histogram::compute_threshold (int min_auto_threshold)
{
	//printf("compute threshold \n");
	u_int64_t sum_allk = 0 ;

	if (_length >= 2)
	{
		_histogram_smoothed [1 ].abundance = (u_int64_t) (  0.6 * (double)  _histogram[1].abundance + 0.4 * (double) _histogram[2].abundance) ;
		sum_allk += _histogram[1].abundance * 1 ;

	}
	
	int index_first_increase = -1;
	int index_maxval_after_first_increase = -1;
	u_int64_t max_val = 0 ;

	//smoothing and detection of first increase
	for (size_t i=2; i<_length  ; i++) // && i < 100
	{
		//printf("idx %i     %i \n",_histogram[i].index ,_histogram[i].abundance );

		sum_allk += _histogram[i].abundance * i  ;
		_histogram_smoothed [i].abundance = // _histogram[i].abundance ;
		(u_int64_t) (
		0.2 * (double) _histogram[i-1].abundance
		+0.6 * (double) _histogram[i].abundance
												+ 0.2 * (double) _histogram[i+1].abundance);
	
		if( index_first_increase==-1 && (_histogram_smoothed[i-1].abundance < _histogram_smoothed[i].abundance) )
		{
			index_first_increase = i-1;
		}
		if(index_first_increase>0 &&  (_histogram_smoothed[i].abundance > max_val))
		{
			max_val = _histogram_smoothed[i].abundance ;
			index_maxval_after_first_increase= i;
		}
	}
	
	sum_allk += _histogram[_length].abundance  *  _length ;

	if(index_first_increase ==-1 )
	{
		_cutoff = min_auto_threshold; //def val
		return;
	}
	
	_firstPeak = index_maxval_after_first_increase;

	DEBUG (("index first increase %i  idx maxval %i \n",index_first_increase,index_maxval_after_first_increase));
	
	u_int64_t min_val = 10000000000LL;
	
	int index_minval = -1;
	
	for (int i=index_first_increase; i<= index_maxval_after_first_increase   ; i++)
	{
		if(_histogram_smoothed[i].abundance < min_val)
		{
			min_val = _histogram_smoothed[i].abundance;
			index_minval = i;
		}
	}
	
	if(index_minval !=-1)
		_cutoff = index_minval;

	u_int64_t sum_elim = 0 ;
	double ratio = 0.0;
	int max_cutoff=0;
	for (size_t i=0; i<_length+1; i++)
	{
		sum_elim +=  _histogram[i].abundance * i ;
		ratio = (double)sum_elim / sum_allk; // ratio elim for cutoff i+1

		DEBUG (("thre %i : %lli elim / %lli   : ratio %f \n",i,sum_elim,sum_allk,ratio ));

		if(ratio >= 0.25)
		{
			max_cutoff = i+1;
			break;
		}
	}
	
	if (_cutoff > max_cutoff)
		_cutoff = max_cutoff;
	
	if (_cutoff< min_auto_threshold)
		_cutoff = min_auto_threshold;

	DEBUG (("cutoff  %i  maxcutoff %i \n",index_minval,max_cutoff));

	/*
	printf("raw values \n");
	for (size_t i=1; i<_length  && i < 100 ; i++)
	{
		printf("idx %i     %i \n",_histogram[i].index ,_histogram[i].abundance );
	}
	
	printf("smoothed values \n");

	for (size_t i=1; i<_length  && i < 100 ; i++)
	{
		printf("idx %i     %i \n",_histogram_smoothed[i].index ,_histogram_smoothed[i].abundance );
		
	}
	*/
	
	_nbsolids =0;
	for (size_t i=_cutoff; i<_length+1   ; i++)
	{
		_nbsolids += _histogram[i].abundance;
		
	}
	
	u_int64_t vol_weak=0;
	u_int64_t volume_total=0;
	
	for (size_t i=0; i<_cutoff ; i++)
	{
		vol_weak += _histogram[i].abundance *i;
	}
	
	for (size_t i=0; i<_length+1   ; i++)
	{
		volume_total += _histogram[i].abundance *i;
	}
	_ratio_weak_volume = (float)vol_weak / (float)volume_total;

}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/
