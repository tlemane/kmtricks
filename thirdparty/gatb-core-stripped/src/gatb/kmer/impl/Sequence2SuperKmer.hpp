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

/** \file Sequence2SuperKmer.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief
 */

#ifndef _SEQUENCE_2_SUPERKMER_HPP_
#define _SEQUENCE_2_SUPERKMER_HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/bank/api/IBank.hpp>
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

/* This functor class takes a Sequence as input and split it into super kmers.
 * Each time a new superkmer is found, the 'processSuperkmer' method is called.
 *
 * NOTE : 'processSuperkmer' is pure virtual and so Sequence2SuperKmer class has to be
 * inherited for providing actual implementation of 'processSuperkmer'.
 *
 * NOTE : 'processSuperkmer' is virtual and called many times, which implies a
 * time overhead. Preliminary tests seem to show that this overhead is acceptable;
 * otherwise a template based approach could be used instead (ie Sequence2SuperKmer is
 * templated by a functor that does the job done in 'processSuperkmer').
 */
template<size_t span>
class Sequence2SuperKmer
{
public:
    /** Shortcut. */
    typedef typename Kmer<span>::Type                                       Type;
    typedef typename Kmer<span>::ModelDirect                                ModelDirect;
    typedef typename Kmer<span>::ModelCanonical                             ModelCanonical;
#ifdef NONCANONICAL
    typedef typename Kmer<span>::template ModelMinimizer <ModelDirect>   Model;
#else
    typedef typename Kmer<span>::template ModelMinimizer <ModelCanonical>   Model;
#endif
    typedef typename Model::Kmer                                            KmerType;
    typedef typename Kmer<span>::Count                                      Count;
    typedef typename Kmer<span>::SuperKmer                                  SuperKmer;

    static const u_int64_t DEFAULT_MINIMIZER = 1000000000 ;


	template<typename KType>
	struct KmerFunctor
	{
		
		Sequence2SuperKmer * caller;
		SuperKmer& superKmer;
		int maxs;
		KmerFunctor (Sequence2SuperKmer * Seq2SupCaller,SuperKmer& superKmer,int maxs) :
		caller(Seq2SupCaller),superKmer(superKmer), maxs(maxs) {}
		
		void operator() (const KType& kmer, size_t idx)  {
			
		 if (kmer.isValid() == false)
		 {
			// printf("non valid kmer %s \n", kmer.value().toString(31).c_str());
			 //caller->_model.toString(kmer)
			 // on invalid kmer : output previous superk utput prev
			 caller->processSuperkmer (superKmer);
			 superKmer.reset();
			 
			 superKmer.minimizer = DEFAULT_MINIMIZER;  //marking will have to restart 'from new'
			 
			 caller->_bankStatsLocal.kmersNbInvalid ++;
			 
			 return;
		 }
			
			caller->_bankStatsLocal.kmersNbValid ++;
			
			/** We get the value of the current minimizer. */
			u_int64_t h = kmer.minimizer().value().getVal();
			
			if(DEFAULT_MINIMIZER == h)
			{
				printf("__ non valid kmer %s \n", kmer.value().toString(31).c_str());

			}
			
			/** We have to set minimizer value if not defined. */
			if (superKmer.isValid() == false)  {  superKmer.minimizer = h;  }
			
			/** If the current super kmer is finished (or max size reached), we dump it. */
			if (h != superKmer.minimizer || superKmer.size() >= (size_t)maxs)
			{
				caller->processSuperkmer (superKmer);
				superKmer.reset();
			}
			
			/** We update the superkmer properties (minimizer value and kmers range). */
			superKmer.minimizer    = h;
			superKmer.addKmer(kmer);
			
		}
	};
	
	
	
    void operator() (bank::Sequence& sequence)
    {
        /** We update statistics about the bank. */
        _bankStatsLocal.update (sequence);

        /** We first check whether we got kmers from the sequence or not. */
		int32_t nbKmers = sequence.getData().size() - _model.getKmerSize() + 1;
		if (nbKmers <= 0)  { return ; }
		
		int maxs = std::min((int)((Type::getSize() - 8 )/2),255) ;  // 8 is because  8 bit used for size of superkmers, not mini size and 255 : max superk size on 8 bits

        /** We create a superkmer object. */
        SuperKmer superKmer (_kmersize, _miniSize);

		//iteration et traitement au fil de l'eau, without large kmer buffer (only small buffer for superkmer now )
		_model.iterate(sequence.getData(), KmerFunctor<KmerType>(this,superKmer,maxs));
		
        //output last superK
        processSuperkmer (superKmer);

        if (_nbWrittenKmers > 500000)   {  _progress->inc (_nbWrittenKmers);  _nbWrittenKmers = 0;  }
    }

    /** Constructor. */
    Sequence2SuperKmer (
        Model&                       model,
        size_t                       nbPasses,
        size_t                       currentPass,
        size_t                       nbPartitions,
        tools::dp::IteratorListener* progress,
        BankStats&                   bankStats
    )
    : _model(model), _pass(currentPass), _nbPass(nbPasses), _nbPartitions(nbPartitions),
      _progress (progress), _nbWrittenKmers(0), _nbSuperKmers(0),
      _bankStatsGlobal(bankStats)
    {
        /** Shortcuts. */
        _kmersize = model.getKmerSize();
        _miniSize = model.getMmersModel().getKmerSize();
    }

    /** Destructor (virtual). */
    virtual ~Sequence2SuperKmer ()
    {
        /** In case we have several passes, we must update sequence information only for first pass. */
        if (_pass==0)  { _bankStatsGlobal += _bankStatsLocal;  }
    }

protected:

    Model&           _model;
    size_t           _pass;
    size_t           _nbPass;
    size_t           _nbPartitions;
    size_t           _kmersize;
    size_t           _miniSize;
    tools::dp::IteratorListener* _progress;
    size_t           _nbWrittenKmers;
    size_t           _nbSuperKmers;
    BankStats&       _bankStatsGlobal;
    BankStats        _bankStatsLocal;

    /** Primitive of the template method operator() */
    virtual void processSuperkmer (SuperKmer& superKmer) { _nbSuperKmers++; }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _SEQUENCE_2_SUPERKMER_HPP_ */

