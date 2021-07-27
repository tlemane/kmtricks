#include <gatb/kmer/impl/HashSorting.hpp>
#include <gatb/tools/collections/impl/OAHash.hpp>
#include <gatb/tools/collections/impl/Hash16.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>
#include <chrono>
#include <sabuhash/sabuhash.hpp>
using namespace std::chrono;
/////////////////////////
// SuperKToHashCommand //
/////////////////////////

namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {

template<size_t span>
SuperKToHashCommand<span>::SuperKToHashCommand(
  SuperKmerBinFiles*  superKstorage,
  int                 fileId,
  int                 kmerSize,
  uint64_t*           r_idx,
  Type*               hash_array,
  uint64_t            window,
  bool                sabuhash
)
  :  _superKstorage(superKstorage), _fileId(fileId), _buffer(0), _buffer_size(0),
    _kmerSize(kmerSize), _r_idx(r_idx), _array(hash_array), _win_size(window), _sabuhash(sabuhash)

{
  Type un;
  un.setVal(1);
  _kmerMask = (un << (_kmerSize * 2)) - 1;
  //_mask_radix.setVal((int64_t) 255);
  //_mask_radix = _mask_radix << ((_kmerSize - 4) * 2);
  _shift      = 2 * (_kmerSize - 1);
  //_shift_radix = ((_kmerSize - 4) * 2);
}

template<size_t span>
inline uint64_t SuperKToHashCommand<span>::hash(Type k)
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
  return (hash % _win_size) + (_fileId * _win_size);
}

template<size_t span>
void SuperKToHashCommand<span>::execute()
{
  uint nbbreak = 0;
  unsigned int nb_bytes_read;
  SabuHash hasher(_kmerSize);
  uint64_t val;
  while (_superKstorage->readBlock(&_buffer, &_buffer_size, &nb_bytes_read, _fileId))
  {
    unsigned char *ptr = _buffer;
    uint8_t nbK;
    int nbsuperkmer_read = 0;
    uint8_t newbyte = 0;

    while (ptr < (_buffer+nb_bytes_read))
    {
      nbK = *ptr; ptr++;
      int rem_size = _kmerSize;
      Type Tnewbyte;
      int nbr = 0;
      _seedk.setVal(0);
      while(rem_size >= 4)
      {
        newbyte = *ptr; ptr++;
        Tnewbyte.setVal(newbyte);
        _seedk = _seedk | (Tnewbyte << (8*nbr));
        rem_size -= 4;
        nbr++;
      }
      int uid = 4;
      if (rem_size > 0)
      {
        newbyte = *ptr ; ptr++;
        Tnewbyte.setVal(newbyte);
        _seedk = (_seedk | (Tnewbyte << (8*nbr)));
        uid = rem_size;
      }
      _seedk = _seedk & _kmerMask;

      uint8_t rem = nbK;
      Type temp = _seedk;
      Type rev_temp = revcomp(temp, _kmerSize);
      Type mink, newnt;
      Type kinsert;
      uint64_t idx;

      bool which = (temp < rev_temp);
      mink = which ? temp : rev_temp;

      idx = *_r_idx;
      if (_sabuhash)
      {
        val = mink.getVal();
        kinsert.setVal((hasher.hash(&val) % _win_size) + (_win_size*_fileId));
      }
      else
        kinsert.setVal(hash(mink));
      _array[idx] = kinsert;
      ++*_r_idx;

      for (int i=0; i<nbK; i++, rem--)
      {
        nbbreak++;
        //bool which = (temp < rev_temp);
        //mink = which ? temp : rev_temp;

        //kinsert.setVal( (hasher.hash((uint64_t*)&mink) % _win_size) + (_fileId*_win_size));
        //kinsert.setVal(hash(mink));
        
        //idx = __sync_fetch_and_add(_r_idx, 1);
        //idx = *_r_idx;
        //_array[idx] = kinsert;
        //++*_r_idx;

        if (rem < 2) break;
        
        if (uid>=4)
        {
          newbyte = *ptr ; ptr++;
          Tnewbyte.setVal(newbyte);
          uid = 0;
        }

        newnt = (Tnewbyte >> (2*uid)) & 3; 
        uid++;
        temp = ((temp << 2) | newnt) & _kmerMask;
        newnt.setVal(comp_NT[newnt.getVal()]);
        rev_temp = ((rev_temp >> 2) | (newnt << _shift)) & _kmerMask;

        which = (temp < rev_temp);
        mink = which ? temp : rev_temp;

        idx = *_r_idx;
        if (_sabuhash)
        {
          val = mink.getVal();
          kinsert.setVal((hasher.hash(&val) % _win_size) + (_win_size*_fileId));
        }
        else
          kinsert.setVal(hash(mink));
        _array[idx] = kinsert;
        ++*_r_idx;
      }
      
      //bool which = (temp < rev_temp);
      //mink = which ? temp : rev_temp;
      ////idx = __sync_fetch_and_add(_r_idx, 1);
      //idx = *_r_idx;
      //kinsert.setVal(hash(mink));
      //_array[idx] = kinsert;
      //++*_r_idx;


      nbsuperkmer_read++;
    }
  }
  if (_buffer != 0) free(_buffer);
}

template<size_t span>
HashSortCommand<span>::HashSortCommand(Type* hash_array, uint64_t array_size)
  : _array(hash_array), _size(array_size)
{

}

template<size_t span>
void HashSortCommand<span>::execute()
{
  std::sort(&_array[0], &_array[_size]);
}

template <size_t span>
HashSortingCommand<span>::HashSortingCommand(
  CountProcessor*       processor,
  size_t                cacheSize,
  IteratorListener*     progress,
  TimeInfo&             timeInfo,
  PartiInfo<5>&         pInfo,
  int                   passi,
  int                   parti,
  size_t                nbCores,
  size_t                kmerSize,
  MemAllocator&         pool,
  vector<size_t>&       offsets,
  SuperKmerBinFiles*    superKstorage,
  uint64_t              window,
  bool                  sabuhash)
  : PartitionsCommand_kx1<span>(
        processor, cacheSize, progress, timeInfo,
        pInfo, passi, parti, nbCores, kmerSize, pool, superKstorage),
    _r_idx(0), _nbItemsPerBankPerPart(offsets), _win_size(window), _sabuhash(sabuhash)
{
  _dispatcher = new Dispatcher(this->_nbCores);
}

template<size_t span>
HashSortingCommand<span>::~HashSortingCommand() {}

template <size_t span>
void HashSortingCommand<span>::execute()
{
  this->_processor->beginPart(
      this->_pass_num, this->_parti_num, this->_cacheSize, this->getName());

  if (this->_superKstorage->getNbItems(this->_parti_num) == 0)
  {
    return;
  }

  //_r_idx = (uint64_t *)CALLOC(256, sizeof(uint64_t));
  _r_idx = (uint64_t*)CALLOC(1, sizeof(uint64_t));
  
  executeRead ();
  executeSort ();
  executeDump ();
    
  FREE(_r_idx);
//  FREE(_array);
  this->_progress->inc(this->_pInfo.getNbKmer(this->_parti_num));
  this->_processor->endPart(this->_pass_num, this->_parti_num);
}

template <size_t span>
void HashSortingCommand<span>::executeRead()
{
  TIME_INFO(this->_timeInfo, "1.read");
  this->_superKstorage->openFile("r", this->_parti_num);

  uint64_t sum_nbKmer = 0;

  {
    LocalSynchronizer synchro(this->_pool.getSynchro());
    this->_pool.align(16);

    size_t nbKmers = this->_pInfo.getNbKmer(this->_parti_num);
    _array = (Type*)this->_pool.pool_malloc(nbKmers*sizeof(Type));
    //_array = (Type*)calloc(nbKmers, sizeof(Type));
  }

  vector<ICommand *> cmds;
  for (size_t i = 0; i < this->_nbCores; i++)
  {
    cmds.push_back(new SuperKToHashCommand<span>(
        this->_superKstorage,
        this->_parti_num,
        this->_kmerSize,
        _r_idx,
        _array,
        _win_size,
        _sabuhash));
  }
  _dispatcher->dispatchCommands(cmds, 0);

  this->_superKstorage->closeFile(this->_parti_num);
}

template<size_t span>
void HashSortingCommand<span>::executeSort()
{
  vector<ICommand*> cmds;
  uint64_t size = *_r_idx;
  
  cmds.push_back(new HashSortCommand<span> (_array, size));

  _dispatcher->dispatchCommands(cmds, 0);
}

template<size_t span>
void HashSortingCommand<span>::executeDump()
{
  //Type first = _array[0];
  size_t nbBanks = _nbItemsPerBankPerPart.size();
  if (nbBanks == 0) nbBanks = 1;

  CounterBuilder solidCounter(nbBanks);
  Type previous_kmer;
  previous_kmer = _array[0];
  solidCounter.init(0);


  for (uint64_t i=1; i<*_r_idx; i++)
  {
    if (previous_kmer != _array[i])
    {
      this->insert(previous_kmer, solidCounter);
      //solidCounter.increase(0);
      previous_kmer = _array[i];
      solidCounter.init(0);

    }
    else
    {
      solidCounter.increase(0);
    }
  }
  this->insert(previous_kmer, solidCounter);
}

  template<size_t span>
  HashByHashCommand<span>:: HashByHashCommand (
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
    tools::storage::impl::SuperKmerBinFiles* 		superKstorage,
    uint64_t window_size,
    bool     sabuhash

  )
    : PartitionsCommand<span> (/*partition,*/ processor, cacheSize, progress, timeInfo, pInfo, passi, parti,nbCores,kmerSize,pool,superKstorage),
      _hashMemory(hashMemory), _window_size(window_size), _fileId(parti), _sabuhash(sabuhash)
  {
  }

  template<size_t span>
  uint64_t HashByHashCommand<span>::hash(Type k)
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

    //template<size_t span>
    //class TempCountFileMerger
    //{
    //  typedef typename Kmer<span>::Type  Type;
    //  typedef tools::misc::Abundance<Type> abundance_t;
    //  typedef std::pair< int , Type> ptcf; //  id pointer , kmer value
    //  struct ptcfcomp { bool operator() (ptcf l,ptcf r) { return ((r.second) < (l.second)); } } ;

    //public:
    //  TempCountFileMerger(int reduceTarget, int chunksize) :_reduceTarget(reduceTarget), _chunksize(chunksize),_idx(0)
    //  {
    //  }

    //  std::vector<string>  mergeFiles(std::vector<string> filenames)
    //  {
    //    ptcf best_elem;
    //    int best_p;
    //    int current_ab = 0;
    //    int previous_ab = 0;
    //    Type current_kmer,previous_kmer;


    //    while(filenames.size() > _reduceTarget)
    //    {

    //      std::vector<string> currentFiles;
    //      for(int ii=0; ii<_chunksize; ii++)
    //      {
    //        currentFiles.push_back(filenames.back()); filenames.pop_back();
    //      }

    //      //the new file containing the merged counts
    //      std::string newfname = currentFiles[0]  + Stringify::format ("_merged_%i", _idx++) ;
    //      BagFile<abundance_t> * bagf = new BagFile<abundance_t>(newfname); LOCAL(bagf);
    //      Bag<abundance_t> * currentbag =  new BagCache<abundance_t> (  bagf, 10000 ); LOCAL(currentbag);

    //      filenames.push_back(newfname);

    //      std::vector<Iterator<abundance_t>*> _tmpCountIterators;

    //      for(int ii=0; ii< currentFiles.size(); ii++)
    //      {
    //        _tmpCountIterators.push_back( new IteratorFile<abundance_t> (currentFiles[ii])  );
    //      }
    //      std::priority_queue< ptcf, std::vector<ptcf>,ptcfcomp > pq;


    //      //// init all iterators  ////
    //      for(int ii=0; ii< _tmpCountIterators.size(); ii++)
    //      {
    //        _tmpCountIterators[ii]->first();
    //      }

    //      //////   init pq ////
    //      for(int ii=0; ii< _tmpCountIterators.size(); ii++)
    //      {
    //        if( ! _tmpCountIterators[ii]->isDone())  {
    //          pq.push(ptcf(ii,_tmpCountIterators[ii]->item().value) );
    //        }
    //      }

    //      //now merge the n sorted iterators and merge their kmer counts.
    //      if(pq.size() != 0)
    //      {
    //        //get first pointer
    //        best_elem = pq.top() ; pq.pop();
    //        best_p = best_elem.first;
    //        previous_ab = _tmpCountIterators[best_p]->item().abundance;
    //        previous_kmer = best_elem.second;

    //        //go forward in this list
    //        _tmpCountIterators[best_p]->next();
    //        if (! _tmpCountIterators[best_p]->isDone())
    //        {
    //          pq.push(ptcf( best_p,_tmpCountIterators[best_p]->item().value) );
    //        }

    //        while (pq.size() != 0)
    //        {

    //          //get  first pointer
    //          best_elem = pq.top() ; pq.pop();
    //          best_p = best_elem.first;
    //          current_ab = _tmpCountIterators[best_p]->item().abundance;
    //          current_kmer = best_elem.second;

    //          //go forward in this list
    //          _tmpCountIterators[best_p]->next();
    //          if (! _tmpCountIterators[best_p]->isDone())
    //          {
    //            pq.push(ptcf( best_p,_tmpCountIterators[best_p]->item().value) );
    //          }


    //          if(current_kmer != previous_kmer)
    //          {
    //            //output previous kmer
    //            currentbag->insert( abundance_t(previous_kmer,previous_ab) );
    //            previous_kmer = current_kmer;
    //            previous_ab = current_ab;
    //          }
    //          else
    //          {
    //            //merge counter
    //            previous_ab += current_ab;
    //          }
    //        }

    //        //output last one
    //        currentbag->insert( abundance_t(previous_kmer,previous_ab) );
    //      }


    //      currentbag->flush();


    //      //cleanup

    //      for(int ii=0; ii< _tmpCountIterators.size(); ii++)
    //      {
    //        delete _tmpCountIterators[ii];
    //      }


    //      //erase used files
    //      for(int ii=0; ii< currentFiles.size(); ii++)
    //      {
    //        std::string fname = currentFiles[ii];
    //        system::impl::System::file().remove(fname);
    //      }

    //    }


    //    return filenames;
    //  }

    //private :

    //  int _reduceTarget;
    //  int _chunksize;
    //  int _idx;

    //};


/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
    template<size_t span>
    void HashByHashCommand<span>:: execute ()
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


      SabuHash hasher(this->_kmerSize);
      uint64_t val;
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
            if (_sabuhash)
            {
              val = mink.getVal();
              mink.setVal((hasher.hash(&val) % _window_size) + (_window_size*_fileId));
            }
            else
              mink.setVal(hash(mink));
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



} } } } // end of namespaces
