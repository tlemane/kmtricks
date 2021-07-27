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

#include <gatb/kmer/impl/PartitionsCommand.hpp>
#include <gatb/tools/collections/impl/OAHash.hpp>
#include <gatb/tools/collections/impl/Hash16.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>
#include <chrono>




using namespace std;
using namespace std::chrono;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::kmer::impl;

#define DEBUG(a)  // printf a

/********************************************************************************/
namespace gatb  {  namespace core  {   namespace kmer  {   namespace impl {
/********************************************************************************/

#define IX(x,rad) ((rad)+(256)*(x))

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
PartitionsCommand<span>:: PartitionsCommand (
   // Iterable<Type>&     partition,
    CountProcessor*     processor,
    size_t              cacheSize,
    IteratorListener*   progress,
    TimeInfo&           timeInfo,
    PartiInfo<5>&       pInfo,
    int                 passi,
    int                 parti,
    size_t              nbCores,
    size_t              kmerSize,
    MemAllocator&       pool,
	tools::storage::impl::SuperKmerBinFiles* 		superKstorage

)
    :
     // _partition(partition),
      _progress(progress),
      _pInfo(pInfo),
      _pass_num(passi),
      _parti_num(parti),
      _nbCores(nbCores),
      _kmerSize(kmerSize),
      _cacheSize(cacheSize),
      _pool(pool),
      _globalTimeInfo(timeInfo),
      _processor(0),
	  _superKstorage(superKstorage)
{
    setProcessor      (processor);
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
PartitionsCommand<span>::~PartitionsCommand()
{
    _globalTimeInfo += _timeInfo;

    setProcessor (0);
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
void PartitionsCommand<span>::insert (const Type& kmer, const CounterBuilder& counter)
{
    /** We call the count processor instance with the information collected for the current kmer. */
    _processor->process (_parti_num, kmer, counter.get());
}

template<size_t span>
PartitionsCommand_kx1<span>:: PartitionsCommand_kx1 (
   // Iterable<Type>&     partition,
    CountProcessor*     processor,
    size_t              cacheSize,
    IteratorListener*   progress,
    TimeInfo&           timeInfo,
    PartiInfo<5>&       pInfo,
    int                 passi,
    int                 parti,
    size_t              nbCores,
    size_t              kmerSize,
    MemAllocator&       pool,
	tools::storage::impl::SuperKmerBinFiles* 		superKstorage

)
    :
     // _partition(partition),
      _progress(progress),
      _pInfo(pInfo),
      _pass_num(passi),
      _parti_num(parti),
      _nbCores(nbCores),
      _kmerSize(kmerSize),
      _cacheSize(cacheSize),
      _pool(pool),
      _globalTimeInfo(timeInfo),
      _processor(0),
	  _superKstorage(superKstorage)
{
    setProcessor      (processor);
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
PartitionsCommand_kx1<span>::~PartitionsCommand_kx1()
{
    _globalTimeInfo += _timeInfo;

    setProcessor (0);
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
void PartitionsCommand_kx1<span>::insert (const Type& kmer, const CounterBuilder& counter)
{
    /** We call the count processor instance with the information collected for the current kmer. */
    _processor->process (_parti_num, kmer, counter.get());
}
	

/////////   multibank version  with old partition system//////////////
template<size_t span>
PartitionsCommand_multibank<span>:: PartitionsCommand_multibank (
											 Iterable<Type>&     partition,
											 CountProcessor*     processor,
											 size_t              cacheSize,
											 IteratorListener*   progress,
											 TimeInfo&           timeInfo,
											 PartiInfo<5>&       pInfo,
											 int                 passi,
											 int                 parti,
											 size_t              nbCores,
											 size_t              kmerSize,
											 MemAllocator&       pool
											 )
:
_partition(partition),
_progress(progress),
_pInfo(pInfo),
_pass_num(passi),
_parti_num(parti),
_nbCores(nbCores),
_kmerSize(kmerSize),
_cacheSize(cacheSize),
_pool(pool),
_globalTimeInfo(timeInfo),
_processor(0)
{
	setProcessor      (processor);
}

template<size_t span>
PartitionsCommand_multibank<span>::~PartitionsCommand_multibank()
{
	_globalTimeInfo += _timeInfo;
	
	setProcessor (0);
}

template<size_t span>
void PartitionsCommand_multibank<span>::insert (const Type& kmer, const CounterBuilder& counter)
{
	/** We call the count processor instance with the information collected for the current kmer. */
	_processor->process (_parti_num, kmer, counter.get());
}
/////////////////
	
	
	
/*********************************************************************
                #     #     #      #####   #     #
                #     #    # #    #     #  #     #
                #     #   #   #   #        #     #
                #######  #     #   #####   #######
                #     #  #######        #  #     #
                #     #  #     #  #     #  #     #
                #     #  #     #   #####   #     #
*********************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
/** in this scheme we count k-mers inside a partition by a hash table */
template<size_t span>
PartitionsByHashCommand<span>:: PartitionsByHashCommand (
   // Iterable<Type>&         partition,
    CountProcessor*         processor,
    size_t                  cacheSize,
    IteratorListener*       progress,
    TimeInfo&               timeInfo,
    PartiInfo<5>&           pInfo,
    int                     passi,
    int                     parti,
    size_t                  nbCores,
    size_t                  kmerSize,
    MemAllocator&           pool,
    u_int64_t               hashMemory,
	tools::storage::impl::SuperKmerBinFiles* 		superKstorage
	//uint64_t window_size

)
    : PartitionsCommand<span> (/*partition,*/ processor, cacheSize, progress, timeInfo, pInfo, passi, parti,nbCores,kmerSize,pool,superKstorage),
     _hashMemory(hashMemory), _window_size(1000), _fileId(parti)
{
}

template<size_t span>
uint64_t PartitionsByHashCommand<span>::hash(Type k)
{
  uint64_t hash = 0;
  uint64_t key = k.getVal();
  hash ^= (hash <<  7) ^  key * (hash >> 3) ^ (~((hash << 11) + (key ^ (hash >> 5))));
  hash = (~hash) + (hash << 21); // hash = (hash << 21) - hash - 1;
  hash = hash ^ (hash >> 24);
  hash = (hash + (hash << 3)) + (hash << 8); // hash * 265
  hash = hash ^ (hash >> 14);
  hash = (hash + (hash << 2)) + (hash << 4); // hash * 21
  hash = hash ^ (hash >> 28);
  hash = hash + (hash << 31);
  return (hash % _window_size) + (_fileId * _window_size);
}
	
//will take N sorted files, will merge them to M files, by chunks of T files at a time


template<size_t span>
TempCountFileMerger<span>::TempCountFileMerger(int reduceTarget, int chunksize)
  : _reduceTarget(reduceTarget), _chunksize(chunksize), _idx(0)
{

}

template<size_t span>
std::vector<string> TempCountFileMerger<span>::mergeFiles(std::vector<string> filenames)
{
  ptcf best_elem;
  int best_p;
  int current_ab = 0;
  int previous_ab = 0;
  Type current_kmer,previous_kmer;


  while(filenames.size() > _reduceTarget)
  {

    std::vector<string> currentFiles;
    for(int ii=0; ii<_chunksize; ii++)
    {
      currentFiles.push_back(filenames.back()); filenames.pop_back();
    }

    //the new file containing the merged counts
    std::string newfname = currentFiles[0]  + Stringify::format ("_merged_%i", _idx++) ;
    BagFile<abundance_t> * bagf = new BagFile<abundance_t>(newfname); LOCAL(bagf);
    Bag<abundance_t> * currentbag =  new BagCache<abundance_t> (  bagf, 10000 ); LOCAL(currentbag);

    filenames.push_back(newfname);

    std::vector<Iterator<abundance_t>*> _tmpCountIterators;

    for(int ii=0; ii< currentFiles.size(); ii++)
    {
      _tmpCountIterators.push_back( new IteratorFile<abundance_t> (currentFiles[ii])  );
    }
    std::priority_queue< ptcf, std::vector<ptcf>,ptcfcomp > pq;


    //// init all iterators  ////
    for(int ii=0; ii< _tmpCountIterators.size(); ii++)
    {
      _tmpCountIterators[ii]->first();
    }

    //////   init pq ////
    for(int ii=0; ii< _tmpCountIterators.size(); ii++)
    {
      if( ! _tmpCountIterators[ii]->isDone())  {
        pq.push(ptcf(ii,_tmpCountIterators[ii]->item().value) );
      }
    }

    //now merge the n sorted iterators and merge their kmer counts.
    if(pq.size() != 0)
    {
      //get first pointer
      best_elem = pq.top() ; pq.pop();
      best_p = best_elem.first;
      previous_ab = _tmpCountIterators[best_p]->item().abundance;
      previous_kmer = best_elem.second;

      //go forward in this list
      _tmpCountIterators[best_p]->next();
      if (! _tmpCountIterators[best_p]->isDone())
      {
        pq.push(ptcf( best_p,_tmpCountIterators[best_p]->item().value) );
      }

      while (pq.size() != 0)
      {

        //get  first pointer
        best_elem = pq.top() ; pq.pop();
        best_p = best_elem.first;
        current_ab = _tmpCountIterators[best_p]->item().abundance;
        current_kmer = best_elem.second;

        //go forward in this list
        _tmpCountIterators[best_p]->next();
        if (! _tmpCountIterators[best_p]->isDone())
        {
          pq.push(ptcf( best_p,_tmpCountIterators[best_p]->item().value) );
        }


        if(current_kmer != previous_kmer)
        {
          //output previous kmer
          currentbag->insert( abundance_t(previous_kmer,previous_ab) );
          previous_kmer = current_kmer;
          previous_ab = current_ab;
        }
        else
        {
          //merge counter
          previous_ab += current_ab;
        }
      }

      //output last one
      currentbag->insert( abundance_t(previous_kmer,previous_ab) );
    }


    currentbag->flush();


    //cleanup

    for(int ii=0; ii< _tmpCountIterators.size(); ii++)
    {
      delete _tmpCountIterators[ii];
    }


    //erase used files
    for(int ii=0; ii< currentFiles.size(); ii++)
    {
      std::string fname = currentFiles[ii];
      system::impl::System::file().remove(fname);
    }

  }


  return filenames;
}


//	template<size_t span>
//	class TempCountFileMerger
//	{
//		typedef typename Kmer<span>::Type  Type;
//		typedef tools::misc::Abundance<Type> abundance_t;
//		typedef std::pair< int , Type> ptcf; //  id pointer , kmer value
//		struct ptcfcomp { bool operator() (ptcf l,ptcf r) { return ((r.second) < (l.second)); } } ;
//
//	public:
//		TempCountFileMerger(int reduceTarget, int chunksize) :_reduceTarget(reduceTarget), _chunksize(chunksize),_idx(0)
//		{
//		}
//
//		std::vector<string>  mergeFiles(std::vector<string> filenames)
//		{
//			ptcf best_elem;
//			int best_p;
//			int current_ab = 0;
//			int previous_ab = 0;
//			Type current_kmer,previous_kmer;
//
//
//			while(filenames.size() > _reduceTarget)
//			{
//
//				std::vector<string> currentFiles;
//				for(int ii=0; ii<_chunksize; ii++)
//				{
//					currentFiles.push_back(filenames.back()); filenames.pop_back();
//				}
//
//				//the new file containing the merged counts
//				std::string newfname = currentFiles[0]  + Stringify::format ("_merged_%i", _idx++) ;
//				BagFile<abundance_t> * bagf = new BagFile<abundance_t>(newfname); LOCAL(bagf);
//				Bag<abundance_t> * currentbag =  new BagCache<abundance_t> (  bagf, 10000 ); LOCAL(currentbag);
//
//				filenames.push_back(newfname);
//
//				std::vector<Iterator<abundance_t>*> _tmpCountIterators;
//
//				for(int ii=0; ii< currentFiles.size(); ii++)
//				{
//					_tmpCountIterators.push_back( new IteratorFile<abundance_t> (currentFiles[ii])  );
//				}
//				std::priority_queue< ptcf, std::vector<ptcf>,ptcfcomp > pq;
//
//
//				//// init all iterators  ////
//				for(int ii=0; ii< _tmpCountIterators.size(); ii++)
//				{
//					_tmpCountIterators[ii]->first();
//				}
//
//				//////   init pq ////
//				for(int ii=0; ii< _tmpCountIterators.size(); ii++)
//				{
//					if( ! _tmpCountIterators[ii]->isDone())  {
//						pq.push(ptcf(ii,_tmpCountIterators[ii]->item().value) );
//					}
//				}
//
//				//now merge the n sorted iterators and merge their kmer counts.
//				if(pq.size() != 0)
//				{
//					//get first pointer
//					best_elem = pq.top() ; pq.pop();
//					best_p = best_elem.first;
//					previous_ab = _tmpCountIterators[best_p]->item().abundance;
//					previous_kmer = best_elem.second;
//
//					//go forward in this list
//					_tmpCountIterators[best_p]->next();
//					if (! _tmpCountIterators[best_p]->isDone())
//					{
//						pq.push(ptcf( best_p,_tmpCountIterators[best_p]->item().value) );
//					}
//
//					while (pq.size() != 0)
//					{
//
//						//get  first pointer
//						best_elem = pq.top() ; pq.pop();
//						best_p = best_elem.first;
//						current_ab = _tmpCountIterators[best_p]->item().abundance;
//						current_kmer = best_elem.second;
//
//						//go forward in this list
//						_tmpCountIterators[best_p]->next();
//						if (! _tmpCountIterators[best_p]->isDone())
//						{
//							pq.push(ptcf( best_p,_tmpCountIterators[best_p]->item().value) );
//						}
//
//
//						if(current_kmer != previous_kmer)
//						{
//							//output previous kmer
//							currentbag->insert( abundance_t(previous_kmer,previous_ab) );
//							previous_kmer = current_kmer;
//							previous_ab = current_ab;
//						}
//						else
//						{
//							//merge counter
//							previous_ab += current_ab;
//						}
//					}
//
//					//output last one
//					currentbag->insert( abundance_t(previous_kmer,previous_ab) );
//				}
//
//
//				currentbag->flush();
//
//
//				//cleanup
//
//				for(int ii=0; ii< _tmpCountIterators.size(); ii++)
//				{
//					delete _tmpCountIterators[ii];
//				}
//
//
//				//erase used files
//				for(int ii=0; ii< currentFiles.size(); ii++)
//				{
//					std::string fname = currentFiles[ii];
//					system::impl::System::file().remove(fname);
//				}
//
//			}
//
//
//			return filenames;
//		}
//
//		private :
//
//		int _reduceTarget;
//		int _chunksize;
//		int _idx;
//
//	};
	
	
/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void PartitionsByHashCommand<span>:: execute ()
{
	typedef typename tools::collections::impl::Hash16<Type>::cell cell_t;

		this->_superKstorage->openFile("r",this->_parti_num);

	this->_processor->beginPart (this->_pass_num, this->_parti_num, this->_cacheSize, this->getName());

	CounterBuilder solidCounter;

	/** We need a map for storing part of solid kmers. */
	//OAHash<Type> hash (_hashMemory);

	Hash16<Type> hash16 (_hashMemory/MBYTE); // now use hash 16 to ensure always finish. needs more ram than OAHash but seems faster




	// If the partition holds kmers (and not superkmers), it would be :
	//      for (it->first(); !it->isDone(); it->next())   {  hash.increment (it->item());  }

	DEBUG (("PartitionsByHashCommand::execute:  fillsolid parti num %i  by oahash --- mem %llu  MB\n",
			this->_parti_num,_hashMemory/MBYTE
));


	typedef tools::misc::Abundance<Type> abundance_t;
	std::vector<string> _tmpCountFileNames;


		//with decompactage
		//superk
		int ks = this->_kmerSize;
		Type un; un.setVal(1);
		//size_t _shift_val = Type::getSize() -8;
		Type kmerMask = (un << (ks*2)) - un;
		size_t shift = 2*(ks-1);

		Type _seedk;

		int _fileId = this->_parti_num;
		unsigned char * _buffer = 0 ;
		unsigned int _buffer_size = 0;




		unsigned int nb_bytes_read;
		while(this->_superKstorage->readBlock(&_buffer, &_buffer_size, &nb_bytes_read, _fileId))
		{
			unsigned char * ptr = _buffer;
			u_int8_t nbK; //number of kmers in the superkmer
			u_int8_t newbyte=0;

			while(ptr < (_buffer+nb_bytes_read)) //decode whole block
			{
				//decode a superkmer
				nbK = *ptr; ptr++;
				//int nb_bytes_superk = (this->_kmerSize + nbK -1 +3) /4  ;

				int rem_size = this->_kmerSize;

				Type Tnewbyte;
				int nbr=0;
				_seedk.setVal(0);
				while(rem_size>=4)
				{
					newbyte = *ptr ; ptr++;
					Tnewbyte.setVal(newbyte);
					_seedk =  _seedk  |  (Tnewbyte  << (8*nbr)) ;
					rem_size -= 4; nbr++;
				}

				int uid = 4; //uid = nb nt used in current newbyte

				//reste du seed kmer
				if(rem_size>0)
				{
					newbyte = *ptr ; ptr++;
					Tnewbyte.setVal(newbyte);

					_seedk = ( _seedk  |  (Tnewbyte  << (8*nbr)) ) ;
					uid = rem_size;
				}
				_seedk = _seedk & kmerMask;



				u_int8_t rem = nbK;
				Type temp = _seedk;
				Type rev_temp = revcomp(temp,this->_kmerSize);
				Type newnt ;
				Type mink;
				Type kinsert;

				//iterate over kmers of this superk
				for (int ii=0; ii< nbK; ii++,rem--)
				{

#ifdef NONCANONICAL
					mink = temp;
#else
					mink = std::min (rev_temp, temp);
#endif

					//mink.setVal(hash(mink));
					/** We insert the kmer into the hash. */
					hash16.insert(mink);


					if(rem < 2) break; //no more kmers in this superkmer, the last one has just been eaten

					////////now decode next kmer of this superkmer ///////

					if(uid>=4) //read next byte
					{
						newbyte = *ptr ; ptr++;
						Tnewbyte.setVal(newbyte);
						uid =0;
					}

					newnt = (Tnewbyte >> (2*uid))& 3; uid++;
					temp = ((temp << 2 ) |  newnt   ) & kmerMask;

					newnt.setVal(comp_NT[newnt.getVal()]) ;
					rev_temp = ((rev_temp >> 2 ) |  (newnt << shift) ) & kmerMask;
				}

				//now go to next superk of this block, ptr should point to beginning of next superk
			}

			//check if hashtable is getting too big : in that case dump it disk and resume with the emptied hashtable
			//at the end merge-sort all the dumped files with the content of hash table
			if(hash16.getByteSize() > _hashMemory) // to be improved (can be slightly larger than maxmemory by a block size)
				//if(_tmpCountFileNames.size()<20) //force  dumps for testing
			{
				//printf("splitting into subparts %lli KB / %lli KB  parti %i subpart %i \n",hash16.getByteSize()/1024,_hashMemory/1024 ,this->_parti_num,_tmpCountFiles.size()  );
				//dump partial count to disk file


				Iterator < cell_t >* itKmerAbundancePartial = hash16.iterator(true);
				LOCAL (itKmerAbundancePartial);


				std::string fname = this->_superKstorage->getFileName(this->_parti_num) + Stringify::format ("_subpart_%i", _tmpCountFileNames.size()) ;
				_tmpCountFileNames.push_back(fname);

				BagFile<abundance_t> * bagf = new BagFile<abundance_t>(fname); LOCAL(bagf);
				Bag<abundance_t> * currentbag =  new BagCache<abundance_t> (  bagf, 10000 ); LOCAL(currentbag);


				for (itKmerAbundancePartial->first(); !itKmerAbundancePartial->isDone(); itKmerAbundancePartial->next())
				{
					cell_t & cell = itKmerAbundancePartial->item();
					currentbag->insert( abundance_t(cell.graine,cell.val) );
				}

				currentbag->flush();
				hash16.clear();
			}
		}


		if(_buffer!=0)
			free(_buffer);



	/** We loop over the solid kmers map.
	 * NOTE !!! we want the items to be sorted by kmer values (see finalize part of debloom). */
	//Iterator < Abundance<Type> >* itKmerAbundance = hash.iterator(true);
	//shortcut
	Iterator < cell_t >* itKmerAbundance = hash16.iterator(true);
	LOCAL (itKmerAbundance);


	//now merge sort over current hash and over the sorted _tmpCountFiles
	// : simple merge sort of n sorted iterators

	if(_tmpCountFileNames.size()!=0)
	{

		TempCountFileMerger<span> tempCountFileMerger (10,10);
		//will merge by chunk of 10 files at a time, until reach less than 10 files
		_tmpCountFileNames = tempCountFileMerger.mergeFiles(_tmpCountFileNames);
		//then will use code below to merge remaining files with the contents of the hash table

		std::vector<Iterator<abundance_t>*> _tmpCountIterators;

		//how to make sure there are not too many subpart files ?  and that we'll not reach the max open files limit ?
		//we *could* merge  only some of them at a time ..  todo ?  --> done with TempCountFileMerger above
		for(int ii=0; ii< _tmpCountFileNames.size(); ii++)
		{
			std::string fname = _tmpCountFileNames[ii];
			_tmpCountIterators.push_back( new IteratorFile<abundance_t> (fname)  );
		}

		// Note (guillaume) : code below is ugly because I have to manage itKmerAbundance (iterator over cell_t)
		// and _tmpCountIterators (iterators over abundance_t) differently since they have different types
		// I would have liked to transform  Iterator<cell_t>  to an Iterator<abundance_t>   with the following adaptor :
		//
		// 		struct cell2AbAdaptor  {  abundance_t operator() (cell_t& c)  { return abundance_t(c.graine,c.val) ; }  };
		//Iterator<abundance_t>*   hashAbundance  = new IteratorAdaptor<cell_t,abundance_t,cell2AbAdaptor> (itKmerAbundance);
	    //	but it turns out to be impossible because of the return by reference of the   item()  function.
		//  Another solution would be to dump contents of hash to a file then read it, but inefficient
		//  So, ugly code it is.  (see all the if(best_p==-1) below)


		//setup the priority queue for merge sorting
		typedef std::pair< int , Type> ptcf; //  id pointer , kmer value
		struct ptcfcomp { bool operator() (ptcf l,ptcf r) { return ((r.second) < (l.second)); } } ;
		std::priority_queue< ptcf, std::vector<ptcf>,ptcfcomp > pq;

		//// init all iterators  ////
		itKmerAbundance->first();
		for(int ii=0; ii< _tmpCountIterators.size(); ii++)
		{
			_tmpCountIterators[ii]->first();
		}

		//////   init pq ////

		if(!itKmerAbundance->isDone())
		{
			pq.push(ptcf(-1,itKmerAbundance->item().graine) ); // -1  will mean in the  itKmerAbundance
		}

		for(int ii=0; ii< _tmpCountIterators.size(); ii++)
		{
			if( ! _tmpCountIterators[ii]->isDone())  {
				abundance_t &ab = _tmpCountIterators[ii]->item();
				pq.push(ptcf(ii,ab.value) );
			}
		}

		ptcf best_elem;
		int best_p;
		int current_ab = 0;
		int previous_ab = 0;
		Type current_kmer,previous_kmer;


	    //now merge the n sorted iterators and merge their kmer counts.
		if(pq.size() != 0)
		{
			//get first pointer
			best_elem = pq.top() ; pq.pop();
			best_p = best_elem.first;
			if(best_p==-1)
			{
				previous_ab = itKmerAbundance->item().val;
			}
			else
			{
				previous_ab = _tmpCountIterators[best_p]->item().abundance;
			}

			previous_kmer = best_elem.second;

			//go forward in this list

			if(best_p==-1)
			{
				itKmerAbundance->next();
				if (! itKmerAbundance->isDone())
				{
					pq.push(ptcf(-1,itKmerAbundance->item().graine) );
				}
			}
			else
			{
				_tmpCountIterators[best_p]->next();
				if (! _tmpCountIterators[best_p]->isDone())
				{
					pq.push(ptcf( best_p,_tmpCountIterators[best_p]->item().value) );
				}
			}

			while (pq.size() != 0)
			{

				//get  first pointer
				best_elem = pq.top() ; pq.pop();
				best_p = best_elem.first;

				if(best_p==-1)
				{
					current_ab = itKmerAbundance->item().val;
				}
				else
				{
					current_ab = _tmpCountIterators[best_p]->item().abundance;
				}
				current_kmer = best_elem.second;

				//go forward in this list
				if(best_p==-1)
				{
					itKmerAbundance->next();
					if (! itKmerAbundance->isDone())
					{
						pq.push(ptcf(-1,itKmerAbundance->item().graine) );
					}
				}
				else
				{
					_tmpCountIterators[best_p]->next();
					if (! _tmpCountIterators[best_p]->isDone())
					{
						pq.push(ptcf( best_p,_tmpCountIterators[best_p]->item().value) );
					}
				}

				if(current_kmer != previous_kmer)
				{
					//output previous kmer
					solidCounter.set (previous_ab);
					this->insert (previous_kmer, solidCounter);
					previous_kmer = current_kmer;
					previous_ab = current_ab;
				}
				else
				{
					//merge counter
					previous_ab += current_ab;
				}

			}

			//output last one
			solidCounter.set (previous_ab);
			this->insert (previous_kmer, solidCounter);
		}


		//cleanup
		for(int ii=0; ii< _tmpCountIterators.size(); ii++)
		{
			delete _tmpCountIterators[ii];
		}


		//erase sub files
		for(int ii=0; ii< _tmpCountFileNames.size(); ii++)
		{
			std::string fname = _tmpCountFileNames[ii];
			system::impl::System::file().remove(fname);
		}

	}
	else // in that case no merging needed, just iterate the hash table and output kmer counts
	{
		for (itKmerAbundance->first(); !itKmerAbundance->isDone(); itKmerAbundance->next())
		{

			cell_t & cell = itKmerAbundance->item();
			solidCounter.set (cell.val);
			this->insert (cell.graine, solidCounter);
		}
	}


	this->_superKstorage->closeFile(this->_parti_num);

	this->_progress->inc (this->_pInfo.getNbKmer(this->_parti_num) ); // this->_pInfo->getNbKmer(this->_parti_num)  kmers.size()

	this->_processor->endPart (this->_pass_num, this->_parti_num);
};

/*********************************************************************
        #     #  #######   #####   #######  #######  ######
        #     #  #        #     #     #     #     #  #     #
        #     #  #        #           #     #     #  #     #
        #     #  #####    #           #     #     #  ######
         #   #   #        #           #     #     #  #   #
          # #    #        #     #     #     #     #  #    #
           #     #######   #####      #     #######  #     #
*********************************************************************/

template<size_t span>
class SuperKReader
{
	typedef typename Kmer<span>::Type  Type;
public:

	void operator() (Type& elem)
	{
		//reading elems by pairs
		if(_first)
		{
			_superk = elem;
			_first = false;
		}
		else
		{
			_seedk = elem;
			
			Type compactedK;
			
			compactedK =  _superk;
			u_int8_t nbK = (compactedK >> _shift_val).getVal()  & 255; // 8 bits poids fort = cpt
			u_int8_t rem = nbK;
			
			Type temp = _seedk;
			Type rev_temp = revcomp(temp,_kmerSize);
			Type newnt ;
			Type mink, prev_mink; prev_mink.setVal(0);
			uint64_t idx;
			
#ifdef NONCANONICAL
            bool prev_which = true;
#else
			bool prev_which =  (temp < rev_temp );
#endif
            
			int kx_size = -1; //next loop start at ii=0, first kmer will put it at 0
			Type radix_kxmer_forward =  (temp & _mask_radix) >> ((_kmerSize - 4)*2);
			Type  first_revk, kinsert,radix_kxmer;
            first_revk.setVal(0);
			
			if(!prev_which) first_revk = rev_temp;
			
			u_int8_t rid;

			for (int ii=0; ii< nbK; ii++,rem--)
			{
#ifdef NONCANONICAL
                bool which = true;
                mink = temp;
#else
				bool which =  (temp < rev_temp );
				mink = which ? temp : rev_temp;
#endif

				if (which != prev_which || kx_size >= _kx) // kxmer_size = 1
				{
					//output kxmer size kx_size,radix_kxmer
					//kx mer is composed of superKp[ii-1] superKp[ii-2] .. superKp[ii-n] with nb elems  n  == kxmer_size +1  (un seul kmer ==k+0)
					
					if(prev_which)
					{
						radix_kxmer = radix_kxmer_forward;
						kinsert = prev_mink;
					}
					else // si revcomp, le radix du kxmer est le debut du dernier kmer
					{
						//previous mink
						radix_kxmer =  (prev_mink & _mask_radix) >> _shift_radix;
						kinsert = first_revk;
					}
					
					//record kxmer
					rid = radix_kxmer.getVal();
					//idx = _r_idx[IX(kx_size,rid)]++;
					idx = __sync_fetch_and_add( _r_idx +  IX(kx_size,rid) ,1); // si le sync fetch est couteux, faire un mini buffer par thread
					
					_radix_kmers [IX(kx_size,rid)][ idx] = kinsert << ((4-kx_size)*2);  //[kx_size][rid]
					if (_bankIdMatrix)  { _bankIdMatrix[IX(kx_size,rid)][ idx] = _bankId; }
					
					radix_kxmer_forward =  (mink & _mask_radix) >> _shift_radix;
					kx_size =0;
					
					if(!which) first_revk = rev_temp;
				}
				else
				{
					kx_size++;
				}
			
				prev_which = which ;
				prev_mink = mink;
				
				if(rem < 2) break; //no more kmers in this superkmer, the last one has just been eaten
				newnt =  ( _superk >> ( 2*(rem-2)) ) & 3 ;
				
				temp = ((temp << 2 ) |  newnt   ) & _kmerMask;
				newnt.setVal(comp_NT[newnt.getVal()]) ;
				rev_temp = ((rev_temp >> 2 ) |  (newnt << _shift) ) & _kmerMask;
			}
			
			//record last kxmer prev_mink et monk ?
			if(prev_which)
			{
				radix_kxmer = radix_kxmer_forward;
				kinsert = prev_mink;
			}
			else // si revcomp, le radix du kxmer est le debut du dernier kmer
			{
				//previous mink
				radix_kxmer =  (prev_mink & _mask_radix) >> _shift_radix;
				kinsert = first_revk;
			}
			
			//record kxmer
			rid = radix_kxmer.getVal();
			//idx = _r_idx[IX(kx_size,rid)]++;
			idx = __sync_fetch_and_add( _r_idx +  IX(kx_size,rid) ,1); // si le sync fetch est couteux, faire un mini buffer par thread
					

            
//                        if (idx >= _radix_sizes[IX(kx_size,rid)] ||  _radix_sizes[IX(kx_size,rid)]  == 0)
//                        {  cout << "error, accessing _radix_kmers beyond bound" << endl; exit(1); }
//                    cout << "tbl ref " << IX(kx_size,rid) << " idx " << idx << " radix sizes " <<   _radix_sizes[IX(kx_size,rid)]   << endl;
            
                        _radix_kmers [IX(kx_size,rid)][ idx] = kinsert << ((4-kx_size)*2);   // [kx_size][rid]
                    //cout << "went okay " << idx << endl;

			if (_bankIdMatrix)  { _bankIdMatrix[IX(kx_size,rid)][ idx] = _bankId; }

			_first = true;
		}
	}
	
	SuperKReader (size_t kmerSize,  uint64_t * r_idx, Type** radix_kmers, uint64_t* radix_sizes, bank::BankIdType** bankIdMatrix, size_t bankId=0)
	: _kmerSize (kmerSize), _kx(4), _radix_kmers(radix_kmers), _radix_sizes(radix_sizes), _bankIdMatrix(bankIdMatrix), _r_idx (r_idx), _first(true), _bankId(bankId)
	 {
		 Type un;
         un.setVal(1);
		 _kmerMask    = (un << (kmerSize*2)) - 1;
		 _mask_radix.setVal((int64_t) 255);
		 _mask_radix  = _mask_radix << ((_kmerSize - 4)*2);
		 _shift       = 2*(kmerSize-1);
		 _shift_val   = un.getSize() -8;
		 _shift_radix = ((kmerSize - 4)*2); // radix is 4 nt long
	}
	
private :

	size_t _kmerSize;
	size_t _shift ;
	size_t _shift_val ;
	size_t _shift_radix ;
	int    _kx;
	Type** _radix_kmers;
      uint64_t*          _radix_sizes;

	bank::BankIdType** _bankIdMatrix;
	uint64_t* _r_idx ;
	bool _first;
	Type _superk, _seedk;
	Type _radix, _mask_radix ;
	Type _kmerMask;
	size_t _bankId;
};

	
	

//pour l'instant marchera que en mode comptage simple, pas en multi jeu separes
//car le jeu separe necessite de passer un cpt de separation des banques, a verifier, et il faut que les buffer conservent l'ordre d'entree, a verifier aussi
//readcommand pour lecture parallele des parti superkmers
template<size_t span>
class ReadSuperKCommand : public gatb::core::tools::dp::ICommand, public system::SmartPointer
{
	typedef typename Kmer<span>::Type  Type;
	
public:
	ReadSuperKCommand(tools::storage::impl::SuperKmerBinFiles* superKstorage, int fileId, int kmerSize,
					  uint64_t * r_idx, Type** radix_kmers, uint64_t* radix_sizes, bank::BankIdType** bankIdMatrix)
	: _superKstorage(superKstorage), _fileId(fileId),_buffer(0),_buffer_size(0), _kmerSize(kmerSize),_radix_kmers(radix_kmers), _radix_sizes(radix_sizes), _bankIdMatrix(bankIdMatrix), _r_idx (r_idx)
	{
		_kx=4;
		Type un;
		un.setVal(1);
		_kmerMask    = (un << (_kmerSize*2)) - 1;
		_mask_radix.setVal((int64_t) 255);
		_mask_radix  = _mask_radix << ((_kmerSize - 4)*2);
		_shift       = 2*(_kmerSize-1);
		_shift_val   = un.getSize() -8;
		_shift_radix = ((_kmerSize - 4)*2); // radix is 4 nt long
	}
	
	void execute ()
	{
	  uint nbbreak = 0;
		unsigned int nb_bytes_read;
		while(_superKstorage->readBlock(&_buffer, &_buffer_size, &nb_bytes_read, _fileId))
		{
			//decode block and iterate through its superkmers
			unsigned char * ptr = _buffer;
			u_int8_t nbK; //number of kmers in the superkmer
			int nbsuperkmer_read = 0;
			u_int8_t newbyte = 0;
			
			while(ptr < (_buffer+nb_bytes_read)) //decode whole block
			{
				//decode a superkmer
				nbK = *ptr; ptr++;
				//int nb_bytes_superk = (_kmerSize + nbK -1 +3) /4  ;
				
				int rem_size = _kmerSize;
				
				Type Tnewbyte;
				int nbr=0;
				_seedk.setVal(0);
				while(rem_size>=4)
				{
					newbyte = *ptr ; ptr++;
					Tnewbyte.setVal(newbyte);
					
					
					_seedk =  _seedk  |  (Tnewbyte  << (8*nbr)) ;
					rem_size -= 4; nbr++;
				}
				
				int uid = 4; //uid = nb nt used in current newbyte
				
				//reste du seed kmer
				if(rem_size>0)
				{
					newbyte = *ptr ; ptr++;
					Tnewbyte.setVal(newbyte);
					
					_seedk = ( _seedk  |  (Tnewbyte  << (8*nbr)) ) ;
					uid = rem_size;
				}
				_seedk = _seedk & _kmerMask;
				
				
				//std::string pt = _seedk.toString(_kmerSize);
				//printf("seedk \n%s \n",pt.c_str());
				
				///////////////////////// seedk should be ready here , now parse kx-mers ////////////////////////////////
				
				u_int8_t rem = nbK;
				Type temp = _seedk;
				Type rev_temp = revcomp(temp,_kmerSize);
				Type newnt ;
				Type mink, prev_mink; prev_mink.setVal(0);
				uint64_t idx;
				
#ifdef NONCANONICAL
				bool prev_which = true;
#else
				bool prev_which =  (temp < rev_temp );
#endif
				
				int kx_size = -1; //next loop start at ii=0, first kmer will put it at 0
				Type radix_kxmer_forward =  (temp & _mask_radix) >> ((_kmerSize - 4)*2);
				Type  first_revk, kinsert,radix_kxmer;
				first_revk.setVal(0);
				
				if(!prev_which) first_revk = rev_temp;
				
				u_int8_t rid;
				
				for (int ii=0; ii< nbK; ii++,rem--)
				{
#ifdef NONCANONICAL
					bool which = true;
					mink = temp;
#else
					bool which =  (temp < rev_temp );
					mink = which ? temp : rev_temp;
#endif
					nbbreak++;
					if (which != prev_which || kx_size >= _kx) // kxmer_size = 1
					{
						//output kxmer size kx_size,radix_kxmer
						//kx mer is composed of superKp[ii-1] superKp[ii-2] .. superKp[ii-n] with nb elems  n  == kxmer_size +1  (un seul kmer ==k+0)
						
						if(prev_which)
						{
							radix_kxmer = radix_kxmer_forward;
							kinsert = prev_mink;
						}
						else // si revcomp, le radix du kxmer est le debut du dernier kmer
						{
							//previous mink
							radix_kxmer =  (prev_mink & _mask_radix) >> _shift_radix;
							kinsert = first_revk;
						}
						
						//record kxmer
						rid = radix_kxmer.getVal();
						//idx = _r_idx[IX(kx_size,rid)]++;
						idx = __sync_fetch_and_add( _r_idx +  IX(kx_size,rid) ,1); // si le sync fetch est couteux, faire un mini buffer par thread
						
						_radix_kmers [IX(kx_size,rid)][ idx] = kinsert << ((4-kx_size)*2);  //[kx_size][rid]
						if (_bankIdMatrix)  { _bankIdMatrix[IX(kx_size,rid)][ idx] = _bankId; }
						
						radix_kxmer_forward =  (mink & _mask_radix) >> _shift_radix;
						kx_size =0;
						
						if(!which) first_revk = rev_temp;
					}
					else
					{
						kx_size++;
					}
					
					prev_which = which ;
					prev_mink = mink;
					
					if(rem < 2) break; //no more kmers in this superkmer, the last one has just been eaten
					
					//////////////////////////////now decode next kmer of this superkmer //////////////////////////////////////////////
					
					if(uid>=4) //read next byte
					{
						newbyte = *ptr ; ptr++;
						Tnewbyte.setVal(newbyte);
						uid =0;
					}
					
					newnt = (Tnewbyte >> (2*uid))& 3; uid++;
					
					temp = ((temp << 2 ) |  newnt   ) & _kmerMask;
					
					newnt.setVal(comp_NT[newnt.getVal()]) ;
					rev_temp = ((rev_temp >> 2 ) |  (newnt << _shift) ) & _kmerMask;
					
					///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				}
				
				//record last kxmer prev_mink et monk ?
				if(prev_which)
				{
					radix_kxmer = radix_kxmer_forward;
					kinsert = prev_mink;
				}
				else // si revcomp, le radix du kxmer est le debut du dernier kmer
				{
					//previous mink
					radix_kxmer =  (prev_mink & _mask_radix) >> _shift_radix;
					kinsert = first_revk;
				}
				
				//record kxmer
				rid = radix_kxmer.getVal();
				//idx = _r_idx[IX(kx_size,rid)]++;
				idx = __sync_fetch_and_add( _r_idx +  IX(kx_size,rid) ,1); // si le sync fetch est couteux, faire un mini buffer par thread
				
				
				_radix_kmers [IX(kx_size,rid)][ idx] = kinsert << ((4-kx_size)*2);   // [kx_size][rid]
				//cout << "went okay " << idx << endl;
				
				if (_bankIdMatrix)  { _bankIdMatrix[IX(kx_size,rid)][ idx] = _bankId; }
				
				
				
				
				//////////////////////////////////////////////////////////
				//ptr+=nb_bytes_superk;
				
				//now go to next superk of this block, ptr should point to beginning of next superk
				nbsuperkmer_read++;
				/////////
			}
			
			//printf("nb superk in this block %i parti %i\n",nbsuperkmer_read,_fileId);
			
			
		}
		
		if(_buffer!=0)
			free(_buffer);
	}
private:
	tools::storage::impl::SuperKmerBinFiles* _superKstorage;
	int _fileId;
	unsigned char * _buffer;
	unsigned int _buffer_size;
	int _kmerSize;
	int    _kx;
	
	Type** _radix_kmers;
	uint64_t*          _radix_sizes;
	bank::BankIdType** _bankIdMatrix;
	uint64_t* _r_idx ;
	
	Type _superk, _seedk;
	Type _radix, _mask_radix ;
	Type _kmerMask;
	size_t _shift ;
	size_t _shift_val ;
	size_t _shift_radix ;
	size_t _bankId;
	
};
	
/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
/** in this scheme we count k-mers in a partition by sorting a vector*/
template<size_t span>
PartitionsByVectorCommand<span>:: PartitionsByVectorCommand (
  //  Iterable<Type>&     partition,
    CountProcessor*     processor,
    size_t              cacheSize,
    IteratorListener*   progress,
    TimeInfo&           timeInfo,
    PartiInfo<5>&       pInfo,
    int                 passi,
    int                 parti,
    size_t              nbCores,
    size_t              kmerSize,
    MemAllocator&       pool,
    vector<size_t>&     offsets,
	tools::storage::impl::SuperKmerBinFiles* 		superKstorage
)
    : PartitionsCommand<span> (/*partition,*/ processor, cacheSize,  progress, timeInfo, pInfo, passi, parti,nbCores,kmerSize,pool,superKstorage),
        _radix_kmers (0), _bankIdMatrix(0), _radix_sizes(0), _r_idx(0), _nbItemsPerBankPerPart(offsets)
{
    _dispatcher = new Dispatcher (this->_nbCores);
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
PartitionsByVectorCommand<span>:: ~PartitionsByVectorCommand ()
{
    if (_dispatcher)  { delete _dispatcher; }
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
void PartitionsByVectorCommand<span>::execute ()
{
    this->_processor->beginPart (this->_pass_num, this->_parti_num, this->_cacheSize, this->getName());

    /** We check that we got something. */

	if (this->_superKstorage->getNbItems(this->_parti_num) == 0)  {  return;  }


    /** We configure tables. */
    _radix_kmers  = (Type**)     MALLOC (256*(KX+1)*sizeof(Type*)); //make the first dims static ?  5*256
    _radix_sizes  = (uint64_t*)  MALLOC (256*(KX+1)*sizeof(uint64_t));
    _r_idx        = (uint64_t*)  CALLOC (256*(KX+1),sizeof(uint64_t));

    /** We need extra information for kmers counting in case of several input banks. */
    if (_nbItemsPerBankPerPart.size() > 1) { _bankIdMatrix = (bank::BankIdType**) MALLOC (256*(KX+1)*sizeof(bank::BankIdType*)); }
    else                                   { _bankIdMatrix = 0; }

    /** We have 3 phases here: read, sort and dump. */
		
		executeRead ();
		executeSort ();
		executeDump ();

    /** We cleanup tables. */
    FREE (_radix_sizes) ;
    FREE (_radix_kmers);
    FREE (_r_idx);
    if (_bankIdMatrix)  { FREE (_bankIdMatrix); }

    /** We update the progress bar. */
    this->_progress->inc (this->_pInfo.getNbKmer(this->_parti_num) );

    this->_processor->endPart (this->_pass_num, this->_parti_num);
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
void PartitionsByVectorCommand<span>::executeRead ()
{
    TIME_INFO (this->_timeInfo, "1.read");

	this->_superKstorage->openFile("r",this->_parti_num);

    /** Recall that the attribute _offsets has a size equals to the number of banks + 1 as input
     * and for each bank, it holds the number of items found for the currently processed partition.
     *
     *               bank0   bank1   ...   bankI
     *   offsets :    xxx     xxx           xxx
     *               <------------------------->
     *                current partition content
     */

	 DEBUG (("_offsets.size=%d  OFFSETS: ", _nbItemsPerBankPerPart.size() ));
	 for (size_t j=0; j<_nbItemsPerBankPerPart.size(); j++)  {  DEBUG (("%6d ", _nbItemsPerBankPerPart[j]));  }  DEBUG (("\n"));

    uint64_t sum_nbxmer =0;

    /** We synchronize this statements block because of threads concurrent access. */
    {
        LocalSynchronizer synchro (this->_pool.getSynchro());

        /** We align the pool with a big alignment constraint (see below macos issue with uint128) */
        this->_pool.align (16);

        /** FIRST: allocation for the kmers. */
        for (size_t xx=0; xx< (KX+1); xx++)
        {
            for (int ii=0; ii< 256; ii++)
            {
                /** Shortcut. */
                size_t nbKmers = this->_pInfo.getNbKmer(this->_parti_num,ii,xx);

                //use memory pool here to avoid memory fragmentation
                _radix_kmers  [IX(xx,ii)] = (Type*)     this->_pool.pool_malloc (nbKmers * sizeof(Type),     "kmers alloc");
                _radix_sizes  [IX(xx,ii)] = nbKmers;

                sum_nbxmer +=  nbKmers;
            }
        }

        /** SECOND: allocation for the bank ids if needed.
         * => NEED TO BE DONE AFTER THE KMERS BECAUSE OF MEMORY ALIGNMENT CONCERNS.
         * On MacOs, we got some crashes with uint128 that were not aligned on 16 bytes
         */
        if (_bankIdMatrix)
		{
			throw Exception ("PartitionsByVectorCommand: multi-bank unsupported with new superk storage");

			/*
			for (size_t xx=0; xx< (KX+1); xx++)
			{
				for (int ii=0; ii< 256; ii++)
				{
					size_t nbKmers = this->_pInfo.getNbKmer(this->_parti_num,ii,xx);
					_bankIdMatrix [IX(xx,ii)] = (bank::BankIdType*) this->_pool.pool_malloc (nbKmers * sizeof(bank::BankIdType), "bank ids alloc");
				}
			}
			 */
		}
    }

    DEBUG (("PartitionsByVectorCommand<span>::executeRead:  fillsolid parti num %i  by vector  nb kxmer / nbkmers      %lli / %lli     %f   with %zu nbcores \n",
        this->_parti_num, sum_nbxmer, this->_pInfo.getNbKmer(this->_parti_num),
        (double) sum_nbxmer /  this->_pInfo.getNbKmer(this->_parti_num),this->_nbCores
    ));

    /** HOW TO COUNT KMERS BY SET OF READS ?
     * Now, we are going to read the temporary partition built during the previous phase and fill
     * the _radix_kmers attribute. We also need to know in _radix_kmers what is the contribution of
     * each bank. We therefore need to iterate the current partition by bank (using the information
     * of _offsets). */

    if (_bankIdMatrix)
	{
		throw Exception ("PartitionsByVectorCommand: multi-bank unsupported with new superk storage");
//		/** We create an iterator over all the items. */
//		Iterator<Type>* itGlobal = this->_partition.iterator();
//		LOCAL (itGlobal);
//		
//		/** We iterate the banks. */
//		for (size_t b=0; b<_nbItemsPerBankPerPart.size(); b++)
//		{
//			/** We truncate the global iterator.
//			 * NB : we initialize (ie call 'first') the global iterator only at first call (for b==0). */
//			Iterator<Type>* itLocal = new TruncateIterator<Type> (*itGlobal, _nbItemsPerBankPerPart[b], b==0 ? true : false);
//			LOCAL (itLocal);
//			
//			/** We iterate this local iterator. */
//			_dispatcher->iterate (itLocal, SuperKReader<span>  (this->_kmerSize, _r_idx, _radix_kmers, _radix_sizes, _bankIdMatrix, b), 10000); //must be even , reading by pairs
//		}
//		
//		/** We check that the global iterator is finished. */
//		if (itGlobal->isDone() == false)  { throw Exception ("PartitionsByVectorCommand: iteration should be finished"); }
		
	}
    else
    {
        /** We iterate the superkmers. */

			vector<ICommand*> cmds;
			for (size_t tid=0; tid < this->_nbCores; tid++)
			{
				cmds.push_back(new ReadSuperKCommand<span> (
															this->_superKstorage,
															this->_parti_num,
															this->_kmerSize,
															_r_idx, _radix_kmers, _radix_sizes, 0
															)
							   );
			}
			
			_dispatcher->dispatchCommands (cmds, 0);

		
		

//		printf("-----done ReadSuperKCommand ---\n");

    }
	
		this->_superKstorage->closeFile(this->_parti_num);

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
class SortCommand : public gatb::core::tools::dp::ICommand, public system::SmartPointer
{
public:
    typedef typename Kmer<span>::Type  Type;

    /** Constructor. */
    SortCommand (Type** kmervec, bank::BankIdType** bankIdMatrix, int begin, int end, uint64_t* radix_sizes)
        : _deb(begin), _fin(end), _radix_kmers(kmervec), _bankIdMatrix(bankIdMatrix), _radix_sizes(radix_sizes) {}

    /** */
    void execute ()
    {
        vector<size_t> idx;
        vector<Tmp>    tmp;

        for (int ii=_deb; ii <=_fin; ii++)
        {
            if (_radix_sizes[ii] > 0)
            {
                /** Shortcuts. */
                Type* kmers = _radix_kmers  [ii];

                if (_bankIdMatrix)
                {
                    /** NOT OPTIMAL AT ALL... in particular we have to use 'idx' and 'tmp' vectors
                     * which may use (a lot of ?) memory. */

                    /** Shortcut. */
                    bank::BankIdType* banksId = _bankIdMatrix [ii];

                    /** NOTE: we sort the indexes, not the items. */
                    idx.resize (_radix_sizes[ii]);
                    for (size_t i=0; i<idx.size(); i++)  { idx[i]=i; }

                    std::sort (idx.begin(), idx.end(), Cmp(kmers));

                    /** Now, we have to reorder the two provided vectors with the same order. */
                    tmp.resize (idx.size());
                    for (size_t i=0; i<idx.size(); i++)
                    {
                        tmp[i].kmer = kmers  [idx[i]];
                        tmp[i].id   = banksId[idx[i]];
                    }
                    for (size_t i=0; i<idx.size(); i++)
                    {
                        kmers  [i] = tmp[i].kmer;
                        banksId[i] = tmp[i].id;
                    }
                }
                else
                {
                    std::sort (&kmers[0] , &kmers[ _radix_sizes[ii]]);
                }
            }
        }
    }

private :

    struct Tmp { Type kmer;  bank::BankIdType id;};

    struct Cmp
    {
        Type* _kmers;
        Cmp (Type* kmers) : _kmers(kmers) {}
        bool operator() (size_t a, size_t b)  { return _kmers[a] < _kmers[b]; }
    };

    int        _deb;
    int        _fin;
    Type**     _radix_kmers;
    bank::BankIdType** _bankIdMatrix;
    uint64_t*  _radix_sizes;
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
void PartitionsByVectorCommand<span>::executeSort ()
{
    TIME_INFO (this->_timeInfo, "2.sort");

    vector<ICommand*> cmds;

    int nwork = 256 / this->_nbCores;

    for (size_t xx=0; xx < (KX+1); xx++)
    {
        cmds.clear();

        //fill cmd work vector
        for (size_t tid=0; tid < this->_nbCores; tid++)
        {
            int deb = 0 + tid * nwork;
            int fin = (tid+1) * nwork -1; // thread will do inclusive range [begin -- end ]
            if(tid == this->_nbCores-1)  { fin = 255; }

            // mettre dans le  SortCommand le master radix_kmers et range a traiter
            cmds.push_back (new SortCommand<span> (
                _radix_kmers+ IX(xx,0),
                (_bankIdMatrix ? _bankIdMatrix+ IX(xx,0) : 0),
                deb, fin,
                _radix_sizes + IX(xx,0)
            ));
        }

        _dispatcher->dispatchCommands (cmds, 0);
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
class KxmerPointer
{
public:
    typedef typename Kmer<span>::Type  Type;

    //on lui passe le vector dun kxmer //std::vector<vector<Type> >  & kmervec
    KxmerPointer (
        Type**      kmervec,
        int         prefix_size,
        int         x_size,
        int         min_radix,
        int         max_radix,
        int         kmerSize,
        uint64_t*   radix_sizes,
        bank::BankIdType**  bankIdMatrix
    )
        : _kxmers(kmervec), _bankIdMatrix(0), _radix_sizes(radix_sizes), _cur_idx(-1),
          _low_radix(min_radix),_high_radix(max_radix),
          _prefix_size(prefix_size), _kmerSize(kmerSize),  _x_size(x_size)
    {
        _idx_radix = min_radix;
        Type un; un.setVal( 1);
        _kmerMask = (un << (_kmerSize*2)) - un;

        _shift_size = ( (4 - _prefix_size) *2) ;
        _radixMask.setVal(_idx_radix) ;
        _radixMask = _radixMask << ((_kmerSize-4)*2);
        _radixMask = _radixMask  << (2*_prefix_size)  ;

        if (bankIdMatrix) { _bankIdMatrix = bankIdMatrix + IX(x_size,0); }
    }

    /** */
    inline bool next ()
    {
        _cur_idx++;

        // go to next non empty radix
        while(_idx_radix<= _high_radix && (uint64_t)_cur_idx >=   _radix_sizes[_idx_radix])
        {
            _idx_radix++;
            _cur_idx = 0;
            //update radix mask does not happen often
            _radixMask.setVal(_idx_radix) ;
            _radixMask = _radixMask << ((_kmerSize-4)*2);
            _radixMask = _radixMask  << (2*_prefix_size)  ;
        }

        return (_idx_radix <= _high_radix);
    }

    /** */
    inline Type    value   () const  {  return ( ((_kxmers[_idx_radix][_cur_idx]) >> _shift_size)  |  _radixMask  ) & _kmerMask ;  }

    /** */
    inline bank::BankIdType getBankId () const  {  return _bankIdMatrix ? _bankIdMatrix [_idx_radix][_cur_idx] : 0;  }

private :

    Type**      _kxmers;
    bank::BankIdType**  _bankIdMatrix;
    uint64_t*   _radix_sizes;
    int64_t     _cur_idx;
    Type        _cur_radix;
    Type        _kmerMask;
    Type        _radixMask;
    int         _idx_radix;
    int         _low_radix;
    int         _high_radix;
    int         _shift_size;
    int         _prefix_size;
    int         _kmerSize;
    int         _x_size; //x size of the _kxmersarray
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
void PartitionsByVectorCommand<span>::executeDump ()
{
    TIME_INFO (this->_timeInfo, "3.dump");

    int nbkxpointers = 453; //6 for k1 mer, 27 for k2mer, 112 for k3mer  453 for k4mer
    vector< KxmerPointer<span>*> vec_pointer (nbkxpointers);
    int best_p;

    std::priority_queue< kxp, std::vector<kxp>,kxpcomp > pq;

    size_t nbBanks = _nbItemsPerBankPerPart.size();
    if (nbBanks == 0) nbBanks = 1;

    CounterBuilder solidCounter (nbBanks);

    Type previous_kmer ;

    //init the pointers to the 6 arrays
    int pidx =0;

    ////////////////////////////////////////////////
    ////-------------k0 pointers-----------/////////
    ////////////////////////////////////////////////

    vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(0,0) ,0,0,0,255,this->_kmerSize, _radix_sizes + IX(0,0), _bankIdMatrix); // vec, prefix size, kxsize , radix min, radix max ,ksize

    ////////////////////////////////////////////////
    ////-------------k1 pointers-----------/////////
    ////////////////////////////////////////////////

    //prefix0
    vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(1,0) ,0,1,0,255,this->_kmerSize, _radix_sizes + IX(1, 0), _bankIdMatrix);
    int lowr = 0;
    int maxr = 63;

    //prefix1
    for(unsigned int ii=0; ii<4; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(1,0) ,1,1,lowr,maxr,this->_kmerSize, _radix_sizes + IX(1, 0), _bankIdMatrix);
        lowr += 64;
        maxr += 64;
    }

    ////////////////////////////////////////////////
    ////-------------k2 pointers-----------/////////
    ////////////////////////////////////////////////

    //prefix0
    vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(2,0),0,2,0,255,this->_kmerSize, _radix_sizes + IX(2, 0), _bankIdMatrix);

    //prefix1
    lowr = 0; maxr = 63;
    for(unsigned int ii=0; ii<4; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(2,0),1,2,lowr,maxr,this->_kmerSize, _radix_sizes + IX(2, 0), _bankIdMatrix);
        lowr += 64;
        maxr += 64;
    }

    //prefix2
    lowr = 0; maxr = 15;
    for(unsigned int ii=0; ii<16; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(2,0),2,2,lowr,maxr,this->_kmerSize, _radix_sizes + IX(2, 0), _bankIdMatrix);
        lowr += 16;
        maxr += 16;
    }

    ////////////////////////////////////////////////
    ////-------------k3 pointers-----------/////////
    ////////////////////////////////////////////////

    //prefix0
    vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(3,0),0,3,0,255,this->_kmerSize, _radix_sizes + IX(3, 0), _bankIdMatrix);

    //prefix1
    lowr = 0; maxr = 63;
    for(unsigned int ii=0; ii<4; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(3,0),1,3,lowr,maxr,this->_kmerSize, _radix_sizes + IX(3, 0), _bankIdMatrix);
        lowr += 64;
        maxr += 64;
    }

    //prefix2
    lowr = 0; maxr = 15;
    for(unsigned int ii=0; ii<16; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(3,0),2,3,lowr,maxr,this->_kmerSize, _radix_sizes + IX(3, 0), _bankIdMatrix);
        lowr += 16;
        maxr += 16;
    }

    //prefix3
    lowr = 0; maxr = 3;
    for(unsigned int ii=0; ii<64; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(3,0),3,3,lowr,maxr,this->_kmerSize, _radix_sizes + IX(3, 0), _bankIdMatrix);
        lowr += 4;
        maxr += 4;
    }

    ////////////////////////////////////////////////
    ////-------------k4 pointers-----------/////////
    ////////////////////////////////////////////////

    //prefix0
    vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(4,0),0,4,0,255,this->_kmerSize, _radix_sizes + IX(4, 0), _bankIdMatrix);

    //prefix1
    lowr = 0; maxr = 63;
    for(unsigned int ii=0; ii<4; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(4,0),1,4,lowr,maxr,this->_kmerSize, _radix_sizes + IX(4, 0), _bankIdMatrix);
        lowr += 64;
        maxr += 64;
    }

    //prefix2
    lowr = 0; maxr = 15;
    for(unsigned int ii=0; ii<16; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(4,0),2,4,lowr,maxr,this->_kmerSize, _radix_sizes + IX(4, 0), _bankIdMatrix);
        lowr += 16;
        maxr += 16;
    }

    //prefix3
    lowr = 0; maxr = 3;
    for(unsigned int ii=0; ii<64; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(4,0),3,4,lowr,maxr,this->_kmerSize, _radix_sizes + IX(4, 0), _bankIdMatrix);
        lowr += 4;
        maxr += 4;
    }

    //prefix4
    lowr = 0; maxr = 0;
    for(unsigned int ii=0; ii<256; ii++)
    {
        vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(4,0),4,4,lowr,maxr,this->_kmerSize, _radix_sizes + IX(4, 0), _bankIdMatrix);
        lowr += 1;
        maxr += 1;
    }

    //fill the  priority queue with the first elems
    for (int ii=0; ii<nbkxpointers; ii++)
    {
        if(vec_pointer[ii]->next())  {  pq.push(kxp(ii,vec_pointer[ii]->value()));  }
    }

    if (pq.size() != 0) // everything empty, no kmer at all
    {
        //get first pointer
        best_p = pq.top().first ; pq.pop();

        previous_kmer = vec_pointer[best_p]->value();

        solidCounter.init (vec_pointer[best_p]->getBankId());

        //merge-scan all 'virtual' arrays and output counts
        while (1)
        {
            //go forward in this array or in new array if reaches end of this one
            if (! vec_pointer[best_p]->next())
            {
                //reaches end of one array
                if(pq.size() == 0) break; //everything done

                //otherwise get new best
                best_p = pq.top().first ; pq.pop();
            }

            if (vec_pointer[best_p]->value() != previous_kmer )
            {
                //if diff, changes to new array, get new min pointer
                pq.push(kxp(best_p,vec_pointer[best_p]->value())); //push new val of this pointer in pq, will be counted later

                best_p = pq.top().first ; pq.pop();

                //if new best is diff, this is the end of this kmer
                if(vec_pointer[best_p]->value()!=previous_kmer )
                {
                    this->insert (previous_kmer, solidCounter);

                    solidCounter.init (vec_pointer[best_p]->getBankId());
                    previous_kmer = vec_pointer[best_p]->value();
                }
                else
                {
                    solidCounter.increase (vec_pointer[best_p]->getBankId());
                }
            }
            else
            {
                solidCounter.increase (vec_pointer[best_p]->getBankId());
            }
        }

        //last elem
        this->insert (previous_kmer, solidCounter);
    }

    /** Cleanup. */
    for (int ii=0; ii<nbkxpointers; ii++)  {  delete vec_pointer[ii];  }
}

	
	
//////////////////Multi bank version of partitionbyvectorcommand //////////////////
	
template<size_t span>
PartitionsByVectorCommand_multibank<span>:: PartitionsByVectorCommand_multibank (
																				 Iterable<Type>&     partition,
																				 CountProcessor*     processor,
																				 size_t              cacheSize,
																				 IteratorListener*   progress,
																				 TimeInfo&           timeInfo,
																				 PartiInfo<5>&       pInfo,
																				 int                 passi,
																				 int                 parti,
																				 size_t              nbCores,
																				 size_t              kmerSize,
																				 MemAllocator&       pool,
																				 vector<size_t>&     offsets
																				 )
: PartitionsCommand_multibank<span> (partition, processor, cacheSize,  progress, timeInfo, pInfo, passi, parti,nbCores,kmerSize,pool),
_radix_kmers (0), _bankIdMatrix(0), _radix_sizes(0), _r_idx(0), _nbItemsPerBankPerPart(offsets)
{
	_dispatcher = new Dispatcher (this->_nbCores);
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
PartitionsByVectorCommand_multibank<span>:: ~PartitionsByVectorCommand_multibank ()
{
	if (_dispatcher)  { delete _dispatcher; }
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
void PartitionsByVectorCommand_multibank<span>::execute ()
{
	this->_processor->beginPart (this->_pass_num, this->_parti_num, this->_cacheSize, this->getName());
	
	/** We check that we got something. */

		if (this->_partition.getNbItems() == 0)  {  return;  }
	
	/** We configure tables. */
	_radix_kmers  = (Type**)     MALLOC (256*(KX+1)*sizeof(Type*)); //make the first dims static ?  5*256
	_radix_sizes  = (uint64_t*)  MALLOC (256*(KX+1)*sizeof(uint64_t));
	_r_idx        = (uint64_t*)  CALLOC (256*(KX+1),sizeof(uint64_t));
	
	/** We need extra information for kmers counting in case of several input banks. */
	if (_nbItemsPerBankPerPart.size() > 1) { _bankIdMatrix = (bank::BankIdType**) MALLOC (256*(KX+1)*sizeof(bank::BankIdType*)); }
	else                                   { _bankIdMatrix = 0; }
	
	/** We have 3 phases here: read, sort and dump. */
	executeRead ();
	executeSort ();
	executeDump ();
	
	/** We cleanup tables. */
	FREE (_radix_sizes) ;
	FREE (_radix_kmers);
	FREE (_r_idx);
	if (_bankIdMatrix)  { FREE (_bankIdMatrix); }
	
	/** We update the progress bar. */
	this->_progress->inc (this->_pInfo.getNbKmer(this->_parti_num) );
	
	this->_processor->endPart (this->_pass_num, this->_parti_num);
};

template<size_t span>
void PartitionsByVectorCommand_multibank<span>::executeRead ()
{
	TIME_INFO (this->_timeInfo, "1.read");
	

	
	/** Recall that the attribute _offsets has a size equals to the number of banks + 1 as input
	 * and for each bank, it holds the number of items found for the currently processed partition.
	 *
	 *               bank0   bank1   ...   bankI
	 *   offsets :    xxx     xxx           xxx
	 *               <------------------------->
	 *                current partition content
	 */
	
 DEBUG (("_offsets.size=%d  OFFSETS: ", _nbItemsPerBankPerPart.size() ));
 for (size_t j=0; j<_nbItemsPerBankPerPart.size(); j++)  {  DEBUG (("%6d ", _nbItemsPerBankPerPart[j]));  }  DEBUG (("\n"));
	
	uint64_t sum_nbxmer =0;
	
	/** We synchronize this statements block because of threads concurrent access. */
	{
		LocalSynchronizer synchro (this->_pool.getSynchro());
		
		/** We align the pool with a big alignment constraint (see below macos issue with uint128) */
		this->_pool.align (16);
		
		/** FIRST: allocation for the kmers. */
		for (size_t xx=0; xx< (KX+1); xx++)
		{
			for (int ii=0; ii< 256; ii++)
			{
				/** Shortcut. */
				size_t nbKmers = this->_pInfo.getNbKmer(this->_parti_num,ii,xx);
				
				//use memory pool here to avoid memory fragmentation
				_radix_kmers  [IX(xx,ii)] = (Type*)     this->_pool.pool_malloc (nbKmers * sizeof(Type),     "kmers alloc");
				_radix_sizes  [IX(xx,ii)] = nbKmers;
				
				sum_nbxmer +=  nbKmers;
			}
		}
		
		/** SECOND: allocation for the bank ids if needed.
		 * => NEED TO BE DONE AFTER THE KMERS BECAUSE OF MEMORY ALIGNMENT CONCERNS.
		 * On MacOs, we got some crashes with uint128 that were not aligned on 16 bytes
		 */
		if (_bankIdMatrix)
		{
			
			for (size_t xx=0; xx< (KX+1); xx++)
			{
				for (int ii=0; ii< 256; ii++)
				{
					size_t nbKmers = this->_pInfo.getNbKmer(this->_parti_num,ii,xx);
					_bankIdMatrix [IX(xx,ii)] = (bank::BankIdType*) this->_pool.pool_malloc (nbKmers * sizeof(bank::BankIdType), "bank ids alloc");
				}
			}
		}
	}
	
	DEBUG (("PartitionsByVectorCommand<span>::executeRead:  fillsolid parti num %i  by vector  nb kxmer / nbkmers      %lli / %lli     %f   with %zu nbcores \n",
			this->_parti_num, sum_nbxmer, this->_pInfo.getNbKmer(this->_parti_num),
			(double) sum_nbxmer /  this->_pInfo.getNbKmer(this->_parti_num),this->_nbCores
));
	
	/** HOW TO COUNT KMERS BY SET OF READS ?
	 * Now, we are going to read the temporary partition built during the previous phase and fill
	 * the _radix_kmers attribute. We also need to know in _radix_kmers what is the contribution of
	 * each bank. We therefore need to iterate the current partition by bank (using the information
	 * of _offsets). */
	
	if (_bankIdMatrix)
	{
		/** We create an iterator over all the items. */
		Iterator<Type>* itGlobal = this->_partition.iterator();
		LOCAL (itGlobal);
		
		/** We iterate the banks. */
		for (size_t b=0; b<_nbItemsPerBankPerPart.size(); b++)
		{
			/** We truncate the global iterator.
			 * NB : we initialize (ie call 'first') the global iterator only at first call (for b==0). */
			Iterator<Type>* itLocal = new TruncateIterator<Type> (*itGlobal, _nbItemsPerBankPerPart[b], b==0 ? true : false);
			LOCAL (itLocal);
			
			/** We iterate this local iterator. */
			_dispatcher->iterate (itLocal, SuperKReader<span>  (this->_kmerSize, _r_idx, _radix_kmers, _radix_sizes, _bankIdMatrix, b), 10000); //must be even , reading by pairs
		}
		
		/** We check that the global iterator is finished. */
		if (itGlobal->isDone() == false)  { throw Exception ("PartitionsByVectorCommand: iteration should be finished"); }
	}
	else
	{
		/** We iterate the superkmers. */

			_dispatcher->iterate (this->_partition.iterator(), SuperKReader<span>  (this->_kmerSize, _r_idx, _radix_kmers, _radix_sizes, 0, 0), 10000); //must be even , reading by pairs
		
		
		
		//		printf("-----done ReadSuperKCommand ---\n");
		
	}
	

	
}

template<size_t span>
void PartitionsByVectorCommand_multibank<span>::executeSort ()
{
	TIME_INFO (this->_timeInfo, "2.sort");
	
	vector<ICommand*> cmds;
	
	int nwork = 256 / this->_nbCores;
	
	for (size_t xx=0; xx < (KX+1); xx++)
	{
		cmds.clear();
		
		//fill cmd work vector
		for (size_t tid=0; tid < this->_nbCores; tid++)
		{
			int deb = 0 + tid * nwork;
			int fin = (tid+1) * nwork -1; // thread will do inclusive range [begin -- end ]
			if(tid == this->_nbCores-1)  { fin = 255; }
			
			// mettre dans le  SortCommand le master radix_kmers et range a traiter
			cmds.push_back (new SortCommand<span> (
			_radix_kmers+ IX(xx,0),
			(_bankIdMatrix ? _bankIdMatrix+ IX(xx,0) : 0),
			deb, fin,
			_radix_sizes + IX(xx,0)
												   ));
		}
		
		_dispatcher->dispatchCommands (cmds, 0);
	}
}

template<size_t span>
void PartitionsByVectorCommand_multibank<span>::executeDump ()
{
	TIME_INFO (this->_timeInfo, "3.dump");
	
	int nbkxpointers = 453; //6 for k1 mer, 27 for k2mer, 112 for k3mer  453 for k4mer
	vector< KxmerPointer<span>*> vec_pointer (nbkxpointers);
	int best_p;
	
	std::priority_queue< kxp, std::vector<kxp>,kxpcomp > pq;
	
	size_t nbBanks = _nbItemsPerBankPerPart.size();
	if (nbBanks == 0) nbBanks = 1;
	
	CounterBuilder solidCounter (nbBanks);
	
	Type previous_kmer ;
	
	//init the pointers to the 6 arrays
	int pidx =0;
	
	////////////////////////////////////////////////
	////-------------k0 pointers-----------/////////
	////////////////////////////////////////////////
	
	vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(0,0) ,0,0,0,255,this->_kmerSize, _radix_sizes + IX(0,0), _bankIdMatrix); // vec, prefix size, kxsize , radix min, radix max ,ksize
	
	////////////////////////////////////////////////
	////-------------k1 pointers-----------/////////
	////////////////////////////////////////////////
	
	//prefix0
	vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(1,0) ,0,1,0,255,this->_kmerSize, _radix_sizes + IX(1, 0), _bankIdMatrix);
	int lowr = 0;
	int maxr = 63;
	
	//prefix1
	for(unsigned int ii=0; ii<4; ii++)
	{
		vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(1,0) ,1,1,lowr,maxr,this->_kmerSize, _radix_sizes + IX(1, 0), _bankIdMatrix);
		lowr += 64;
		maxr += 64;
	}
	
	////////////////////////////////////////////////
	////-------------k2 pointers-----------/////////
	////////////////////////////////////////////////
	
	//prefix0
	vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(2,0),0,2,0,255,this->_kmerSize, _radix_sizes + IX(2, 0), _bankIdMatrix);
	
	//prefix1
	lowr = 0; maxr = 63;
	for(unsigned int ii=0; ii<4; ii++)
	{
		vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(2,0),1,2,lowr,maxr,this->_kmerSize, _radix_sizes + IX(2, 0), _bankIdMatrix);
		lowr += 64;
		maxr += 64;
	}
	
	//prefix2
	lowr = 0; maxr = 15;
	for(unsigned int ii=0; ii<16; ii++)
	{
		vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(2,0),2,2,lowr,maxr,this->_kmerSize, _radix_sizes + IX(2, 0), _bankIdMatrix);
		lowr += 16;
		maxr += 16;
	}
	
	////////////////////////////////////////////////
	////-------------k3 pointers-----------/////////
	////////////////////////////////////////////////
	
	//prefix0
	vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(3,0),0,3,0,255,this->_kmerSize, _radix_sizes + IX(3, 0), _bankIdMatrix);
	
	//prefix1
	lowr = 0; maxr = 63;
	for(unsigned int ii=0; ii<4; ii++)
	{
		vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(3,0),1,3,lowr,maxr,this->_kmerSize, _radix_sizes + IX(3, 0), _bankIdMatrix);
		lowr += 64;
		maxr += 64;
	}
	
	//prefix2
	lowr = 0; maxr = 15;
	for(unsigned int ii=0; ii<16; ii++)
	{
		vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(3,0),2,3,lowr,maxr,this->_kmerSize, _radix_sizes + IX(3, 0), _bankIdMatrix);
		lowr += 16;
		maxr += 16;
	}
	
	//prefix3
	lowr = 0; maxr = 3;
	for(unsigned int ii=0; ii<64; ii++)
	{
		vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(3,0),3,3,lowr,maxr,this->_kmerSize, _radix_sizes + IX(3, 0), _bankIdMatrix);
		lowr += 4;
		maxr += 4;
	}
	
	////////////////////////////////////////////////
	////-------------k4 pointers-----------/////////
	////////////////////////////////////////////////
	
	//prefix0
	vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(4,0),0,4,0,255,this->_kmerSize, _radix_sizes + IX(4, 0), _bankIdMatrix);
	
	//prefix1
	lowr = 0; maxr = 63;
	for(unsigned int ii=0; ii<4; ii++)
	{
		vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(4,0),1,4,lowr,maxr,this->_kmerSize, _radix_sizes + IX(4, 0), _bankIdMatrix);
		lowr += 64;
		maxr += 64;
	}
	
	//prefix2
	lowr = 0; maxr = 15;
	for(unsigned int ii=0; ii<16; ii++)
	{
		vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(4,0),2,4,lowr,maxr,this->_kmerSize, _radix_sizes + IX(4, 0), _bankIdMatrix);
		lowr += 16;
		maxr += 16;
	}
	
	//prefix3
	lowr = 0; maxr = 3;
	for(unsigned int ii=0; ii<64; ii++)
	{
		vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(4,0),3,4,lowr,maxr,this->_kmerSize, _radix_sizes + IX(4, 0), _bankIdMatrix);
		lowr += 4;
		maxr += 4;
	}
	
	//prefix4
	lowr = 0; maxr = 0;
	for(unsigned int ii=0; ii<256; ii++)
	{
		vec_pointer[pidx++] = new KxmerPointer<span> (_radix_kmers+ IX(4,0),4,4,lowr,maxr,this->_kmerSize, _radix_sizes + IX(4, 0), _bankIdMatrix);
		lowr += 1;
		maxr += 1;
	}
	
	//fill the  priority queue with the first elems
	for (int ii=0; ii<nbkxpointers; ii++)
	{
		if(vec_pointer[ii]->next())  {  pq.push(kxp(ii,vec_pointer[ii]->value()));  }
	}
	
	if (pq.size() != 0) // everything empty, no kmer at all
	{
		//get first pointer
		best_p = pq.top().first ; pq.pop();
		
		previous_kmer = vec_pointer[best_p]->value();
		
		solidCounter.init (vec_pointer[best_p]->getBankId());
		
		//merge-scan all 'virtual' arrays and output counts
		while (1)
		{
			//go forward in this array or in new array if reaches end of this one
			if (! vec_pointer[best_p]->next())
			{
				//reaches end of one array
				if(pq.size() == 0) break; //everything done
				
				//otherwise get new best
				best_p = pq.top().first ; pq.pop();
			}
			
			if (vec_pointer[best_p]->value() != previous_kmer )
			{
				//if diff, changes to new array, get new min pointer
				pq.push(kxp(best_p,vec_pointer[best_p]->value())); //push new val of this pointer in pq, will be counted later
				
				best_p = pq.top().first ; pq.pop();
				
				//if new best is diff, this is the end of this kmer
				if(vec_pointer[best_p]->value()!=previous_kmer )
				{
					this->insert (previous_kmer, solidCounter);
					
					solidCounter.init (vec_pointer[best_p]->getBankId());
					previous_kmer = vec_pointer[best_p]->value();
				}
				else
				{
					solidCounter.increase (vec_pointer[best_p]->getBankId());
				}
			}
			else
			{
				solidCounter.increase (vec_pointer[best_p]->getBankId());
			}
		}
		
		//last elem
		this->insert (previous_kmer, solidCounter);
	}
	
	/** Cleanup. */
	for (int ii=0; ii<nbkxpointers; ii++)  {  delete vec_pointer[ii];  }
}
	
	
/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
