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

#ifndef _PARTI_INFO_HPP_
#define _PARTI_INFO_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>
#include <mutex>
#include <queue>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

/** class containing info of each parti : exact number of kmers per parti
 *
 * will be computed by fillparti, then used by fillsolids
 * needed beacause the exact number of kmers per parti can no longer be inferred from partition file size
 * xmer ==1 -->  only k0 mer    (xmer == nb max of kmers in kx mer, not the x   k1 mer == 2 kmers max in k1 mer)
 */
template <size_t xmer>
class PartiInfo
{
    struct mmer_bin_record {
        uint64_t nb_superks, nb_kmers, nb_kxmers;
        mmer_bin_record& operator+=(const mmer_bin_record& other) {
            nb_superks += other.nb_superks;
            nb_kmers += other.nb_kmers;
            nb_kxmers += other.nb_kxmers;
            return *this;
        }
    };

    struct parti_record {
        uint64_t nb_kmers, nb_kxmers;
        uint64_t nbk_per_radix[xmer * 256];

        uint64_t& getNbKmer(int radix, int xx)
        { return nbk_per_radix[xx * 256 + radix]; }

        parti_record& operator+=(const parti_record& other) {
            nb_kmers += other.nb_kmers;
            nb_kxmers += other.nb_kxmers;
            for(int i = 0 ; i < xmer * 256 ; i++) {
                nbk_per_radix[i] += other.nbk_per_radix[i];
            }
            return *this;
        }
    };

public:

    inline void incKmer (int numpart, u_int64_t val=1)
    {
        _parti_records[numpart].nb_kmers+=val;
    }

    inline void incKxmer(int numpart, u_int64_t val=1)
    {
        _parti_records[numpart].nb_kxmers+=val; //now used to store number of kxmers per parti
    }

    //superksize in number of kmer inside this superk
    //increases both superk and kmer count
    inline void incSuperKmer_per_minimBin(int numbin, int superksize, u_int64_t val=1)
    {
		_nb_superk_total+=val;
		_nb_kmer_total += val*superksize;
        _mmer_bin_records[numbin].nb_superks+= val ;
        _mmer_bin_records[numbin].nb_kmers+= val *superksize;
    }

    //kxmer count (regardless of x size), used for ram
    inline void incKxmer_per_minimBin(int numbin, u_int64_t val=1)
    {
        _mmer_bin_records[numbin].nb_kxmers+= val ;
    }

    //numaprt, radix, size of kx mer
    inline void incKmer_and_rad(int numpart, int radix,int x,  u_int64_t val=1) //appele ds vrai loop
    {
        parti_record& record = _parti_records[numpart];
        record.nb_kxmers += val;            // now used to store number of kxmers per parti
        record.nb_kmers += (val * (x+1));  // number of  'real' kmers
        record.getNbKmer(radix, x)+=val; // contains number of kx mer per part per radix per x
    }

    /** Adds count to *this from other */
    PartiInfo &add(const PartiInfo &other) {
        for (int np = 0; np < _nbpart; np++) {
            _parti_records[np] += other._parti_records[np];
        }

        for (int ii = 0; ii < _num_mm_bins; ii++) {
            _mmer_bin_records[ii] += other._mmer_bin_records[ii];
        }

        _nb_superk_total += other._nb_superk_total;
        _nb_kmer_total += other._nb_kmer_total;

        return *this;
    }

    /** Adds count to *this from other while preventing concurrent access to *this */
    PartiInfo &add_sync(const PartiInfo &other) {
        // Locking is way faster than doing individuals "lock add"s for all the fields
        std::lock_guard<std::mutex> guard{_mut};
        return add(other);
    }

    /** */
    inline  u_int64_t getNbKmer(int numpart) const
    {
        return _parti_records[numpart].nb_kmers;
    }

    /** Get nbk in bin radix of parti numpart */
    inline  u_int64_t getNbKmer(int numpart, int radix, int xx) const
    {
        return _parti_records[numpart].getNbKmer(radix, xx);
    }

    /** */
    inline  u_int64_t   getNbSuperKmer(int numpart) const //now used for number of kxmers (indistinctive of xsize)
    {
        return _parti_records[numpart].nb_kxmers;
    }
	
	inline  u_int64_t   getNbSuperKmerTotal() const
	{
		return _nb_superk_total;
	}

	inline  u_int64_t   getNbKmerTotal() const
	{
		return _nb_kmer_total;
	}
    /** */
    inline  u_int64_t   getNbSuperKmer_per_minim(int numbin) const
    {
        return _mmer_bin_records[numbin].nb_superks;
    }

    /** */
    inline  u_int64_t   getNbKmer_per_minim(int numbin) const
    {
        return _mmer_bin_records[numbin].nb_kmers;
    }

    /** */
    inline  u_int64_t   getNbKxmer_per_minim(int numbin) const
    {
        return _mmer_bin_records[numbin].nb_kxmers;
    }

    /** */
    void clear()
    {
        memset (_parti_records,    0, _nbpart      * sizeof(parti_record));
        memset (_mmer_bin_records, 0, _num_mm_bins * sizeof(mmer_bin_record));
    }

    /** */
    void printInfo() const
    {
        printf("------------------\n");
        printf("Nb kmers per parti\n");

        for (int np=0; np<_nbpart; np++)  {  printf("Parti[%i]= %lli\n",np,this->getNbKmer(np));  }

        printf("------------------------\n");
        printf("Nb kxmers per parti\n");

        for (int np=0; np<_nbpart; np++)  {  printf("Parti[%i]= %lli\n",np,this->getNbSuperKmer(np));  }

        //			printf("----------------------------\n");
        //			printf("Nb kmers per parti per radix\n");
        //
        //			for (int np=0; np<_nbpart; np++) {
        //				printf("___ Parti %i ___\n",np);
        //
        //				for (int rad=0; rad<256; rad++) {
        //					printf("%10lli  ",this->getNbKmer(np,rad,0));
        //					if((rad & 7) == 7) printf("\n");
        //				}
        //				printf("\n");
        //
        //			}

        printf("----------------------------\n");
        printf("Nb Super kmers , nb kmers per minim bin\n");

        u_int64_t sumk = 0;
        u_int64_t sumsuperk = 0;

        for (int np=0; np<_num_mm_bins; np++)
        {
			if(this->getNbSuperKmer_per_minim(np)!=0 || this->getNbKmer_per_minim(np)  !=0 )
			   	printf("Bin[%5i]= %lli    %lli\n",np, this->getNbSuperKmer_per_minim(np),this->getNbKmer_per_minim(np)    );

            sumk += this->getNbKmer_per_minim(np);
            sumsuperk +=  this->getNbSuperKmer_per_minim(np);
        }

        printf("total number of kmers %lli  total number of superkmers %lli \n",sumk,sumsuperk);
        printf("Average size of superkmers :  %f \n",(double)sumk / (double)sumsuperk );
    }

    /** Constructor. */
    PartiInfo(int nbpart, int minimsize) : _nbpart(nbpart)
    {
		_nb_superk_total =0;
		_nb_kmer_total =0;
        _parti_records = (parti_record*) CALLOC (_nbpart, sizeof(parti_record));
        _num_mm_bins =   1 << (2*minimsize);
        _mmer_bin_records = (mmer_bin_record*) CALLOC (_num_mm_bins, sizeof(mmer_bin_record));
    }


    /** Constructor (copy). Needed for fillparti class. */
    PartiInfo(const PartiInfo& cr)
    //the copy contr realloc its own  arrays, zero init
    {
        _num_mm_bins = cr._num_mm_bins;
        _nbpart      = cr._nbpart;
		_nb_superk_total = cr._nb_superk_total;
		_nb_kmer_total = cr._nb_kmer_total;
		
        _parti_records = (parti_record*) CALLOC (_nbpart, sizeof(parti_record));
        _mmer_bin_records = (mmer_bin_record*) CALLOC (_num_mm_bins, sizeof(mmer_bin_record));

        //	printf("PartiInfo copy constr %p  _nb_kmers_per_parti %p \n",this,_nb_kmers_per_parti);

    }

    /** Destructor. */
    ~PartiInfo()
    {
        FREE (_parti_records);
        FREE (_mmer_bin_records);
    }


    void saveInfoFile(const std::string& prefix)
    {
        std::ofstream myfile;
        myfile.open (prefix + "/PartiInfoFile");
        myfile << _nbpart << "\n";
        myfile << _num_mm_bins << "\n";
        myfile << _nb_superk_total << "\n";
        myfile << _nb_kmer_total << "\n";
        for (int np = 0; np < _nbpart; np++) {
            myfile <<  _parti_records[np].nb_kmers << "\n";
            myfile <<  _parti_records[np].nb_kxmers << "\n";
             for(int i = 0 ; i < xmer * 256 ; i++) {
                myfile <<  _parti_records[np].nbk_per_radix[i] << "\n";
            }
        }
        for (int ii = 0; ii < _num_mm_bins; ii++) {
            myfile << _mmer_bin_records[ii].nb_superks << "\n";
            myfile << _mmer_bin_records[ii].nb_kmers << "\n";
            myfile << _mmer_bin_records[ii].nb_kxmers << "\n";
        }
        myfile.close();
    }
    
    PartiInfo(const std::string& prefix)
    {
        std::ifstream myfile (prefix+"/PartiInfoFile");
        std::string line;
        if (myfile.is_open())
        {
            getline (myfile,line); _nbpart = std::stoull(line);
            getline (myfile,line); _num_mm_bins = std::stoull(line);
            getline (myfile,line); _nb_superk_total = std::stoull(line);
            getline (myfile,line); _nb_kmer_total = std::stoull(line);

            _parti_records = (parti_record*) CALLOC (_nbpart, sizeof(parti_record));
            _mmer_bin_records = (mmer_bin_record*) CALLOC (_num_mm_bins, sizeof(mmer_bin_record));
            for (int np = 0; np < _nbpart; np++) {
                getline (myfile,line); _parti_records[np].nb_kmers = std::stoull(line);
                getline (myfile,line); _parti_records[np].nb_kxmers = std::stoull(line);
                for(int i = 0 ; i < xmer * 256 ; i++) {
                    getline (myfile,line); _parti_records[np].nbk_per_radix[i] = std::stoull(line);
                }
            }
            for (int ii = 0; ii < _num_mm_bins; ii++) {
                getline (myfile,line); _mmer_bin_records[ii].nb_superks = std::stoull(line);
                getline (myfile,line); _mmer_bin_records[ii].nb_kmers = std::stoull(line);
                getline (myfile,line); _mmer_bin_records[ii].nb_kxmers = std::stoull(line);
            }
            myfile.close();
        }
        else
        {
            std::cout << "error reading PartiInfoFile at " << prefix << std::endl; exit(1);
        }
    }


private:
    std::mutex _mut = {}; // Locking for merging with main PartiInfo
    u_int64_t _nbpart;
    u_int64_t _num_mm_bins;
	u_int64_t _nb_superk_total;
	u_int64_t _nb_kmer_total;
    parti_record*  _parti_records;
    mmer_bin_record*  _mmer_bin_records;
};

/********************************************************************************/

/** Class providing hash values for minimizers.
 *
 * Minimizer (of size m) values are gathered in N partitions, so a function from 4^m values
 * to N values is built. This function is built thanks to a PartiInfo object that is built
 * while scanning a part of the input bank (ie. information about minimizers distribution
 * is built then).
 *
 * The hash value is got through the operator() method.
 *
 * It is possible to save/load an instance of this class; this may be useful for having the
 * same partition as the solid kmers for other purpose (like debloom for instance).
 *
 * Note: this class could be renamed in the future to emphasize it essentially provides a
 * hash function service.
 */
class Repartitor : public system::SmartPointer
{
public:

    /** Hash value type. */
    typedef u_int16_t Value;

    /** Constructor
     * \param[in] nbpart : hash value will be in range [0..nbpart-1]
     * \param[in] minimsize : size of the minimizers. */
    Repartitor (int nbpart=0, int minimsize=0, int nbPass=1)  : _nbpart(nbpart), _nb_minims(1 << (minimsize*2)), _nbPass(nbPass), _freq_order(0)
    {
        if (nbpart <= 0)  { system::Exception("Repartitor: nbpart (%d) should be > 0", nbpart); }
    }

    /** Constructor */
    Repartitor (tools::storage::impl::Group& group)  : _nbpart(0), _nb_minims(0), _nbPass(0), _freq_order(0)   { this->load (group);  }

    /** Destructor */
    ~Repartitor ()  {  if (_freq_order)  { delete[] _freq_order; } }

    /** Compute the hash function for the minimizer.
     * \param[in] pInfo : information about the distribution of the minimizers. */
    void computeDistrib (const PartiInfo<5>& pInfo);
    void justGroupNaive (const PartiInfo<5>& pInfo,  std::vector <std::pair<int,int> > &counts);
    void justGroup      (const PartiInfo<5>& pInfo,  std::vector <std::pair<int,int> > &counts);
    void justGroupLexi  (const PartiInfo<5>& extern_pInfo);

    /** Returns the hash value for the given minimizer value.
     * \param[in] minimizerValue : minimizer value as an integer.
     * \return hash value for the given minimizer. */
    Value operator() (u_int64_t minimizerValue)  { return getRepartTable() [minimizerValue]; }

    /** Load the repartition table from a storage object.
     * \param[in] group : group where the repartition table has to be loaded */
    void load (tools::storage::impl::Group& group);

    /** Save the repartition table into a storage object.
     * \param[in] group : group where the repartition table has to be saved */
    void save (tools::storage::impl::Group& group);

    /** For debug purpose. */
    void printInfo ();

    /** Get the number of passes used to split the input bank. */
    size_t getNbPasses() const { return _nbPass; }

    /** Get a buffer on minimizer frequencies. */
    uint32_t* getMinimizerFrequencies () { return _freq_order; }

    /** Set the minimizer frequencies. */
    void setMinimizerFrequencies (uint32_t* freq) { _freq_order = freq; }

private:

    typedef std::pair<u_int64_t,u_int64_t> ipair; //taille bin, numero bin

    class itriple
    {
        public:
            u_int64_t first;
            u_int64_t second;
            u_int64_t third;

            itriple( u_int64_t first,  u_int64_t second,  u_int64_t third) : first(first), second(second), third(third) {}
            itriple() : first(0), second(0), third(0) {}
    };

    struct compBin {
        bool operator() (ipair l,ipair r) { return l.first > r.first; }
    } comp_bins;

    struct compSpace {
        bool operator() (ipair l,ipair r) { return l.second > r.second; } //for partition, pair is parti number, space ued
    } ;

    struct compSpaceTriple {
        bool operator() (itriple l, itriple r) { return l.second > r.second; } // same as compSpace
    } ;

    typedef std::vector<Value> Table;

    /** Get the repartition table. It is built at first call. */
    const Table& getRepartTable()
    {
		// if (_repart_table.empty()) { throw system::Exception ("Repartitor : table has not been initialized"); }
        return _repart_table;
    }

    u_int16_t          _nbpart;
    u_int64_t          _nb_minims;
    u_int16_t          _nbPass;
    std::vector<Value> _repart_table ;

    uint32_t* _freq_order;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _PARTI_INFO_HPP_ */
