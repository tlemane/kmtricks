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

#include <gatb/kmer/impl/PartiInfo.hpp>
#include <algorithm>

// We use the required packages
using namespace std;

#define DEBUG(a) //printf a
// wanted to debug things separately:
// now there are many different and mutually exclusive debug messages
// DEBUG has more output than debug2
// to see top20 repartitions of bins, uncomment IFDEBUG and DEBUG2
#define IFDEBUG(a)  //a
#define DEBUG2(a) //printf a

/********************************************************************************/
namespace gatb  {  namespace core  {   namespace kmer  {   namespace impl {
/********************************************************************************/

static const u_int32_t MAGIC_NUMBER = 0x12345678;

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Repartitor::computeDistrib (const PartiInfo<5>& extern_pInfo)
{
    /** We allocate a table whose size is the number of possible minimizers. */
    _repart_table.resize (_nb_minims);

    std::vector<ipair> bin_size_vec;
    std::priority_queue< itriple, std::vector<itriple>,compSpaceTriple > pq;

    //sum total bins size
    IFDEBUG(u_int64_t sumsizes =0);
    for (u_int64_t ii=0; ii< _nb_minims; ii++)
    {
        // sumsizes +=   extern_pInfo.getNbSuperKmer_per_minim(ii); // _binsize[ii];
        // bin_size_vec.push_back(ipair( extern_pInfo.getNbSuperKmer_per_minim(ii) ,ii));
        IFDEBUG(sumsizes += extern_pInfo.getNbKxmer_per_minim(ii)); // _binsize[ii];
        bin_size_vec.push_back(ipair( extern_pInfo.getNbKxmer_per_minim(ii) ,ii));
    }

    DEBUG(("Repartitor : mean size per parti should be :  %lli  (total %lli )\n", sumsizes / _nbpart, sumsizes));

    //init space left
    for (int jj = 0; jj < _nbpart; jj++)  {  pq.push (itriple(jj,0,0));  }

    //sort minim bins per size
    std::sort (bin_size_vec.begin (), bin_size_vec.end (), comp_bins);

    DEBUG2(("Repartitor : 20 largest estimated bin sizes \n"));
    for (size_t ii=0; ii<20 &&  ii< bin_size_vec.size(); ii++ )
    {
        DEBUG2 (("binsize [%llu] = %llu \n",bin_size_vec[ii].second,bin_size_vec[ii].first));
    }

    //GC suggestion : put the largest in the emptiest (la plus grosse dans la plus grosse)

    itriple smallest_parti;

    u_int64_t cur_minim = 0;
    while (cur_minim < _nb_minims)
    {
        //get emptiest parti
        smallest_parti = pq.top(); pq.pop();

        //put largest bin in it
        _repart_table[bin_size_vec[cur_minim].second] = smallest_parti.first;

        //update space used in this bin, push it back in the pq
        smallest_parti.second += bin_size_vec[cur_minim].first;
        smallest_parti.third ++; // how many minimizers are in this bin (just for info)

        pq.push (smallest_parti);

        DEBUG (("Repartitor : affected minim %llu to part %llu  space used %llu  (msize %llu) \n",
            bin_size_vec[cur_minim].second,smallest_parti.first,
            smallest_parti.second , bin_size_vec[cur_minim].first
        ));

        cur_minim++;
    }
}

// simple version of the code above in the case where we use frequency-based minimizers, and we just want to group minimizers according to their ordering
void Repartitor::justGroupNaive (const PartiInfo<5>& extern_pInfo, std::vector <std::pair<int,int> > &counts)
{
    /** We allocate a table whose size is the number of possible minimizers. */
    _repart_table.resize (_nb_minims);

    for (u_int64_t ii=0; ii< _nb_minims; ii++)
    {
        _repart_table[ii] = 0; // important to have a consistent repartition for unseen (in the sample approximation) minimizers
    }

    int step = counts.size() / _nbpart;
    
    for (unsigned int i = 0; i < counts.size(); i++)
    {
        _repart_table[counts[i].second] = std::min((int)(i / step), _nbpart-1);
    }
    
}


// much more effective version of the function above, using estimation of number of kmers per bucket
void Repartitor::justGroup (const PartiInfo<5>& extern_pInfo, std::vector <std::pair<int,int> > &counts)
{
    /** We allocate a table whose size is the number of possible minimizers. */
    _repart_table.resize (_nb_minims);

    for (u_int64_t ii=0; ii< _nb_minims; ii++)
    {
        _repart_table[ii] = _nbpart - 1; 
        // important to have a consistent repartition for unseen (in the sample approximation) minimizers
        // and since any minimizer is potentially unseen, for bcalm, we have to put them into the last partition
    }


    IFDEBUG(
            {// same debugging as computeDistrib() 
            std::vector<ipair> bin_size_vec;
            //sum total bins size
            u_int64_t sumsizes =0;
            for (u_int64_t ii=0; ii< _nb_minims; ii++)
            {
            sumsizes += extern_pInfo.getNbKxmer_per_minim(ii);
            bin_size_vec.push_back(ipair( extern_pInfo.getNbKxmer_per_minim(ii) ,ii));
            }
            DEBUG2(("Repartitor(justGroup): mean size per parti should be :  %lli  (total %lli )\n", sumsizes / _nbpart, sumsizes));
            //sort minim bins per size
            std::sort (bin_size_vec.begin (), bin_size_vec.end (), comp_bins);
            DEBUG2(("Repartitor(justGroup): 20 largest estimated bin sizes \n"));
            for (size_t ii=0; ii<20 &&  ii< bin_size_vec.size(); ii++ )
            DEBUG2 (("binsize [%llu] = %llu \n",bin_size_vec[ii].second,bin_size_vec[ii].first));
            })

    //sum total count size
    u_int64_t total_counts =0;
    for (unsigned int i = 0; i < counts.size(); i++)
        total_counts += counts[i].first;

    u_int64_t sumsizes =0;
    for (u_int64_t ii=0; ii< _nb_minims; ii++)
        sumsizes += extern_pInfo.getNbKmer_per_minim(ii);
 
    u_int64_t mean_size = sumsizes / _nbpart;
    
    u_int64_t acc = 0, j = 0;
    for (unsigned int i = 0; i < counts.size(); i++)
    {
        _repart_table[counts[i].second] = j;

        acc += extern_pInfo.getNbKmer_per_minim(counts[i].second);
        if (acc > mean_size)
        {
            acc = 0;
            if (j < _nbpart)
                j++;
        }
    }
}

// lexi case
void Repartitor::justGroupLexi (const PartiInfo<5>& extern_pInfo)
{
    /** We allocate a table whose size is the number of possible minimizers. */
    _repart_table.resize (_nb_minims);

    for (u_int64_t ii=0; ii< _nb_minims; ii++)
    {
        _repart_table[ii] = _nbpart-1; 
        // important to have a consistent repartition for unseen (in the sample approximation) minimizers
        // and since any minimizer is potentially unseen, for bcalm, we have to put them into the last partition
    }

    u_int64_t sumsizes =0;
    for (u_int64_t ii=0; ii< _nb_minims; ii++)
        sumsizes += extern_pInfo.getNbKmer_per_minim(ii);
 
    u_int64_t mean_size = sumsizes / _nbpart;
    
    u_int64_t acc = 0, j = 0;
    for (unsigned int i = 0; i < _nb_minims; i++)
    {
        _repart_table[i] = j;
        acc += extern_pInfo.getNbKmer_per_minim(i);
        if (acc > mean_size)
        {
            acc = 0;
            if (j < _nbpart)
                j++;
        }
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
void Repartitor::load (tools::storage::impl::Group& group)
{
    bool hasMinimizerFrequencies = false;

    tools::storage::impl::Storage::istream is (group, "minimRepart");
    is.read ((char*)&_nbpart,     sizeof(_nbpart));
    is.read ((char*)&_nb_minims,  sizeof(_nb_minims));
    is.read ((char*)&_nbPass,     sizeof(_nbPass));

    DEBUG (("[Repartitor::load] :  _nbpart=%d  _nb_minims=%d  _nbPass=%d \n",
        _nbpart, _nb_minims, _nbPass
    ));

    /** We allocate a table whose size is the number of possible minimizers. */
    _repart_table.resize (_nb_minims);

    is.read ((char*)_repart_table.data(), sizeof(Value) * _nb_minims);

    is.read ((char*)&hasMinimizerFrequencies, sizeof(bool));

    u_int32_t magic = 0;
    is.read ((char*)&magic,  sizeof(magic));
    if (magic != MAGIC_NUMBER)  { throw system::Exception("Unable to load Repartitor (minimRepart), possibly due to bad format."); }

    if (hasMinimizerFrequencies)
    {
        tools::storage::impl::Storage::istream is2 (group, "minimFrequency");
        _freq_order = new uint32_t [_nb_minims];
        is2.read ((char*)_freq_order,     sizeof(uint32_t)*_nb_minims);

        is2.read ((char*)&magic,  sizeof(magic));
        if (magic != MAGIC_NUMBER)  { throw system::Exception("Unable to load Repartitor (minimFrequency), possibly due to bad format."); }
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
void Repartitor::save (tools::storage::impl::Group& group)
{
    DEBUG (("[Repartitor::save] :  _nbpart=%d  _nb_minims=%d  _nbPass=%d \n",
        _nbpart, _nb_minims, _nbPass
    ));

    bool hasMinimizerFrequencies = _freq_order != NULL;

    tools::storage::impl::Storage::ostream os (group, "minimRepart");
    os.write ((const char*)&_nbpart,                sizeof(_nbpart));
    os.write ((const char*)&_nb_minims,             sizeof(_nb_minims));
    os.write ((const char*)&_nbPass,                sizeof(_nbPass));
    os.write ((const char*)_repart_table.data(),    sizeof(Value) * _nb_minims);
    os.write ((const char*)&hasMinimizerFrequencies,sizeof(bool));
    os.write ((const char*)&MAGIC_NUMBER,           sizeof(MAGIC_NUMBER));
    os.flush();

    if (hasMinimizerFrequencies)
    {
        tools::storage::impl::Storage::ostream os2 (group, "minimFrequency");
        os2.write ((const char*)_freq_order,    sizeof(uint32_t) * _nb_minims);
        os2.write ((const char*)&MAGIC_NUMBER,  sizeof(MAGIC_NUMBER));
        os2.flush();
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
void Repartitor::printInfo ()
{
    printf("Repartitor : nbMinimizers=%ld\n", _nb_minims);
    for (u_int64_t ii=0; ii<_nb_minims; ii++ )  {  printf("   table[%ld] = %d \n",ii,_repart_table[ii]); }
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
