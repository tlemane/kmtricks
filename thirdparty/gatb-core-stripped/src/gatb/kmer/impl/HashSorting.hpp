#pragma once
#include <gatb/kmer/impl/PartitionsCommand.hpp>
//#include <gatb/kmer/impl/Hashers.hpp>
#include <gatb/kmer/impl/sabuhash.hpp>
//#include <gatb/gatb_core.hpp>
using namespace std;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::collections;
using namespace gatb::core::tools::collections::impl;

using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

using namespace gatb::core::kmer::impl;
using namespace gatb::core::tools::storage::impl;

namespace gatb
{
namespace core
{
namespace kmer
{
namespace impl
{

template <size_t span>
class SuperKToHashCommand : public gatb::core::tools::dp::ICommand, public system::SmartPointer
{
  typedef typename Kmer<span>::Type Type;

public:
  SuperKToHashCommand(
    SuperKmerBinFiles*  superKstorage,
    int                 fileId,
    int                 kmerSize,
    uint64_t*           r_idx,
    Type*               hash_array,
    uint64_t            window,
    bool                sabuhash
  );

  void execute();
  uint64_t hash(Type k);

private: 
  tools::storage::impl::SuperKmerBinFiles *_superKstorage;
  int           _fileId;
  unsigned char *_buffer;
  unsigned int  _buffer_size;
  int           _kmerSize;

  uint64_t  *_r_idx;
  Type      _seedk;
  Type      _kmerMask;
  size_t    _shift;
  Type*     _array;
  uint64_t  _win_size;
  bool      _sabuhash;
};


template<size_t span>
class HashByHashCommand : public PartitionsCommand<span>
{
public:

  /** Shortcut. */ /* R: don't know how to avoid this code duplication => R1: I'm afraid it's not possible. */
  typedef typename Kmer<span>::Type           Type;
  typedef typename Kmer<span>::Count          Count;
  typedef ICountProcessor<span> CountProcessor;

  /** Constructor. */
  HashByHashCommand (
    //    gatb::core::tools::collections::Iterable<Type>& partition,
    CountProcessor*                                 processor,
    size_t                                          cacheSize,
    gatb::core::tools::dp::IteratorListener*        progress,
    tools::misc::impl::TimeInfo&                    timeInfo,
    PartiInfo<5>&                                   pInfo,
    int                                             passi,
    int                                             parti,
    size_t                                          nbCores,
    size_t                                          kmerSize,
    gatb::core::tools::misc::impl::MemAllocator&    pool,
    u_int64_t                                       hashMemory,
    tools::storage::impl::SuperKmerBinFiles* 		superKstorage,
    uint64_t window_size,
    bool      sabuhash
  );

  /** Get the class name (for statistics). */
  const char* getName() const { return "hash"; }

  /** */
  void execute ();
  uint64_t hash(Type k);

private:
  u_int64_t _hashMemory;
  uint64_t _window_size;
  uint64_t _fileId;
  bool     _sabuhash;
};



template <size_t span>
class HashSortingCommand : public PartitionsCommand_kx1<span>
{
public:
  typedef typename Kmer<span>::Type   Type;
  typedef typename Kmer<span>::Count  Count;
  typedef ICountProcessor<span>       CountProcessor;

public:
  HashSortingCommand(
      CountProcessor      *processor,
      size_t              cacheSize,
      IteratorListener    *progress,
      TimeInfo            &timeInfo,
      PartiInfo<5>        &pInfo,
      int                 passi,
      int                 parti,
      size_t              nbCores,
      size_t              kmerSize,
      MemAllocator        &pool,
      vector<size_t>      &offsets,
      SuperKmerBinFiles   *superKstorage,
      uint64_t            window,
      bool                sabuhash);

  ~HashSortingCommand();

  const char *getName() const { return "vector"; }

  void execute();

private:
  uint64_t                *_r_idx;
  uint64_t                _win_size;
  Type                    *_array;
  tools::dp::IDispatcher  *_dispatcher;
  bool                    _sabuhash;

  void executeRead();
  void executeSort();
  void executeDump();

  std::vector<size_t> _nbItemsPerBankPerPart;
};

template<size_t span>
class HashSortCommand : public gatb::core::tools::dp::ICommand, public system::SmartPointer
{
public:
  typedef typename Kmer<span>::Type Type;

  HashSortCommand(Type* hash_array, uint64_t array_size);
  void execute();

private:
  Type* _array;
  uint64_t _size;
};



} } } } // end of namespaces