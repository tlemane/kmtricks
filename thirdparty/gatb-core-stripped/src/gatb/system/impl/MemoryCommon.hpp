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

/** \file MemoryCommon.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Implementations common to various OS.
 */

#ifndef _GATB_CORE_SYSTEM_IMPL_MEMORY_COMMON_HPP_
#define _GATB_CORE_SYSTEM_IMPL_MEMORY_COMMON_HPP_

/********************************************************************************/

#include <gatb/system/api/IMemory.hpp>
#include <gatb/system/api/Exception.hpp>

#include <stdlib.h>
#include <string.h>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace system    {
/** \brief Implementation of Operating System abstraction layer */
namespace impl      {
/********************************************************************************/

/** \brief Implementation of IMemoryAllocator interface using standard system functions.
 *
 *  This implementation provides a few methods common to all operating systems.
 *  It uses functions from stdlib.h
 */
class MemoryAllocatorStdlib : public IMemoryAllocator
{
public:

    /** Singleton. */
    static IMemoryAllocator& singleton()  { static MemoryAllocatorStdlib instance; return instance; }

    /** \copydoc IMemoryAllocator::malloc */
     void* malloc  (BlockSize_t size)
     {
         void* res = ::malloc (size);
         if (!res)  {  throw Exception ("no memory for malloc"); }
         return res;
     }

     /** \copydoc IMemoryAllocator::calloc */
     void* calloc  (size_t nmemb, BlockSize_t size)
     {
         void* res = ::calloc (nmemb, size);
         if (!res)  {  throw Exception ("no memory for calloc"); }
         return res;
     }

     /** \copydoc IMemoryAllocator::realloc */
     void* realloc (void *ptr, BlockSize_t size)
     {
         void* res = ::realloc (ptr, size);
         if (!res)  {  throw Exception ("no memory for realloc"); }
         return res;
     }

     /** \copydoc IMemoryAllocator::free */
     void  free (void *ptr)
     {
         if (ptr != 0)  {  ::free (ptr);  }
     }
};

/********************************************************************************/

/** \brief Implementation of IMemory interface with Proxy design pattern
 *
 * This implementation is a Proxy design pattern: it references a memory allocator
 * and uses this reference for creating memory blocks.
 *
 * Statistics operations are provided through the IMemory interface. Providing such
 * services implies that we have some management of the allocated blocks. In particular
 * we need to store the size of block which is done by allocating some extra space for
 * each block (this space being used for storing the block size).
 *
 * Note that such an implementation requires a little bit more memory than clients
 * actually ask for: for instance, if a client calls malloc with 100, the final allocated
 * block size will be 100+4=104 bytes (if sizeof(BlockSize_t)==4).
 *
 * It is mandatory that blocks created through this implementation are also deleted by the
 * same implementation.
 */
class MemorySizeStore : public IMemoryAllocator
{
public:

    /** Constructor.
     * \param[in] alloc : the referred memory allocator.
     */
    MemorySizeStore (IMemoryAllocator& alloc)
        :  _alloc(alloc), _nbBlocks(0), _currentMemory(0), _peakMemory(0) {}

    /** \copydoc IMemoryAllocator::malloc */
     void* malloc  (BlockSize_t size)
     {
         /** We add the size for storing the size of the pointer to be allocated. */
         BlockSize_t actualSize = size + sizeof(BlockSize_t);

         /** We allocate the block in a classical way. Note the cast due to the fact
          * that we need to do some arithmetic on the pointer. */
         u_int8_t* res = (u_int8_t*) _alloc.malloc (actualSize);

         /** We store the required size. */
         storeBlockSize (res, actualSize);

         /** We update the statistics information. Note that we force synchronization since
          * we could have concurrent access on the instance. */
         __sync_fetch_and_add (&_nbBlocks,      1);
         __sync_fetch_and_add (&_currentMemory, actualSize);

         /** we return the result. */
         return res + sizeof(BlockSize_t);
     }

     /** \copydoc IMemoryAllocator::calloc */
     void* calloc  (size_t nmemb, BlockSize_t size)
     {
         /** We rely on the malloc to emulate the calloc operation. */
         void* res = malloc (nmemb * size);

         /** We reset the allocated memory. */
         ::memset (res, 0, nmemb * size);

         /** We return the result. */
         return res;
     }

     /** \copydoc IMemoryAllocator::realloc */
     void* realloc (void* ptr, BlockSize_t size)
     {
         u_int8_t*   actualPtr    = 0;
         BlockSize_t previousSize = 0;

         /** We retrieve the actual block pointer and the previous block size. */
         if (ptr != 0)
         {
             actualPtr    = (u_int8_t*)ptr - sizeof(BlockSize_t);
             previousSize = *((BlockSize_t*)actualPtr);
         }

         /** We add the size for storing the size of the pointer to be allocated. */
         BlockSize_t actualSize = size + sizeof(BlockSize_t);

         /** We reallocate the actual block. */
         u_int8_t* res = (u_int8_t*)_alloc.realloc (actualPtr, actualSize);

         /** We store the required size. */
         storeBlockSize (res, size);

         /** We update the statistics information. Note that we force synchronization since
          * we could have concurrent access on the instance. */
         if (size > previousSize)  // The size of the block is going to increase.
         {
             __sync_fetch_and_add (&_currentMemory, size - previousSize);
         }
         else  // The size of the block is going to decrease.
         {
             __sync_fetch_and_sub (&_currentMemory, previousSize - size);
         }

         /** We return the result. */
         return res + sizeof(BlockSize_t);
     }

     /** \copydoc IMemoryAllocator::free */
     void  free (void* ptr)
     {
         if (ptr != 0)
         {
             /** We retrieve the actual block pointer. */
             u_int8_t* actualPtr = (u_int8_t*)ptr - sizeof(BlockSize_t);

             /** We retrieve the block size. */
             BlockSize_t actualSize = *((BlockSize_t*)actualPtr);

             /** We release the actual allocated block. */
             _alloc.free (actualPtr);

             /** We update the current blocks number. */
             __sync_fetch_and_sub (&_nbBlocks,      1);
             __sync_fetch_and_sub (&_currentMemory, actualSize);
         }
     }

     /** \copydoc IMemory::getNbBlocks */
     size_t getNbBlocks () { return _nbBlocks; }

     /** \copydoc IMemory::getCurrentUsage */
     TotalSize_t getCurrentUsage () { return _currentMemory; }

     /** \copydoc IMemory::getMaximumUsage */
     TotalSize_t getMaximumUsage () { return _peakMemory; }

protected:

     IMemoryAllocator&  _alloc;

     size_t      _nbBlocks;
     TotalSize_t _currentMemory;
     TotalSize_t _peakMemory;

     /** */
     void storeBlockSize (void* ptr, BlockSize_t size)  {  *((BlockSize_t*)ptr) = size;  }
};

/********************************************************************************/

/** \brief Implementation of IMemory interface with bounded memory usage
 *
 * This implementation is a Proxy design pattern: it references a memory allocator
 * and checks that client requests doesn't exceed this kind of specifications.
 * If the client request is allowed, the request is done by the referred allocator.
 *
 * Two thresholds can be parameterize such an allocator at construction:
 *      - maximum total size allowed
 *      - maximum block size allowed
 *
 * Implementation note: in order to ensure that the thresholds are not exceeded, we need
 * to know the size of each allocated pointer. Therefore, when allocating some block of
 * size N, we actually allocate N+d bytes, where d is the size of an integer big enough
 * to represent the requested block size. This size is provided through the class template T.
 *
 */
class MemoryBounded : public MemorySizeStore
{
public:

    /** Constructor.
     * \param[in] alloc : the referred IMemoryAllocator instance
     * \param[in] maxBlockSize : maximum size allowed for one block allocation
     * \param[in] maxTotalSize : maximum size allowed for the sum of blocks allocation
     * */
    MemoryBounded (IMemoryAllocator& alloc, BlockSize_t maxBlockSize, TotalSize_t maxTotalSize)
    : MemorySizeStore(alloc), _maxBlockSize(maxBlockSize), _maxTotalSize(maxTotalSize) {}

    /** \copydoc IMemoryAllocator::malloc */
     void* malloc  (BlockSize_t size)
     {
         /** We check that the required block size is not too big. */
         if (size > _maxBlockSize)
         {
             throw Exception ("block size too big for malloc: %d required but %d allowed", size, _maxBlockSize);
         }

         /** We check that we don't reach the maximum size allowed. */
         if (_currentMemory + size >= _maxTotalSize)
         {
             throw Exception ("memory maximum reached for malloc: required %d, current %d, max %d", size, _currentMemory, _maxTotalSize);
         }

         /** we return the result. */
         return MemorySizeStore::malloc (size);
     }

     /** \copydoc IMemoryAllocator::calloc */
     void* calloc  (size_t nmemb, BlockSize_t size)
     {
         /** We check that the required block size is not too big. */
         if (nmemb*size > _maxBlockSize)
         {
             throw Exception ("block size too big for calloc: %d required but %d allowed", size, _maxBlockSize);
         }

         /** We check that we don't reach the maximum size allowed. */
         if (_currentMemory + nmemb*size >= _maxTotalSize)
         {
             throw Exception ("memory maximum reached for calloc: required %d, current %d, max %d", size, _currentMemory, _maxTotalSize);
         }

         /** we return the result. */
         return MemorySizeStore::calloc (nmemb, size);
     }

     /** \copydoc IMemoryAllocator::realloc */
     void* realloc (void* ptr, BlockSize_t size)
     {
         /** We check that the required block size is not too big. */
         if (size > _maxBlockSize)
         {
             throw Exception ("block size too big for realloc: %d required but %d allowed", size, _maxBlockSize);
         }

         /** We check that we don't reach the maximum size allowed. */
         if (_currentMemory + size >= _maxTotalSize)
         {
             throw Exception ("memory maximum reached for realloc: required %d, current %d, max %d", size, _currentMemory, _maxTotalSize);
         }

         /** we return the result. */
         return MemorySizeStore::realloc (ptr, size);
     }

     /** \copydoc IMemoryAllocator::free */
     void  free (void* ptr)
     {
         /** We release the actual allocated block. */
         MemorySizeStore::free (ptr);
     }

protected:

     BlockSize_t _maxBlockSize;
     TotalSize_t _maxTotalSize;
};

/********************************************************************************/

/** \brief Implementation of IMemoryOperations interface using standard system functions.
 *
 *  This implementation provides a few methods common to all operating systems.
 *  It uses functions from string.h
 */
class MemoryOperationsCommon : public IMemoryOperations
{
public:

    /** Singleton. */
    static IMemoryOperations& singleton()  { static MemoryOperationsCommon instance; return instance; }

    /** \copydoc IMemoryOperations::memset */
    void* memset (void* s, int c, size_t n)                 { return ::memset (s, c, n);      }

    /** \copydoc IMemoryOperations::memcpy */
    void* memcpy (void* dest, const void* src, size_t n)    { return ::memcpy (dest, src, n); }

    /** \copydoc IMemoryOperations::memcmp */
    int   memcmp (const void* s1, const void* s2, size_t n) { return ::memcmp (s1, s2, n);    }
};

/********************************************************************************/

/** \brief common implementation of IMemory interface
 *
 * This implementation delegates the allocation part to a referred IMemoryAllocator instance.
 *
 * It is not abstract since it implements the statistics methods but with dummy return values.
 * It could have sense in case we want fast allocators without statistic information.
 *
 * Its main purpose is to factorize some code for concrete implementations.
 */
class MemoryCommon : public IMemory
{
public:

    /** Constructor.
     * \param[in] alloc : the referred memory allocator.
     * \param[in] ope   : the referred memory operations.
     */
    MemoryCommon (IMemoryAllocator& alloc, IMemoryOperations& ope) : _alloc(alloc), _ope(ope)  {}

    /** \copydoc IMemoryAllocator::malloc */
     void* malloc  (BlockSize_t size)                { return _alloc.malloc (size);         }

     /** \copydoc IMemoryAllocator::calloc */
     void* calloc  (size_t nmemb, BlockSize_t size)  { return _alloc.calloc (nmemb, size);  }

     /** \copydoc IMemoryAllocator::realloc */
     void* realloc (void *ptr, BlockSize_t size)     { return _alloc.realloc (ptr, size);   }

     /** \copydoc IMemoryAllocator::free */
     void  free (void *ptr)                          { _alloc.free (ptr);                   }

     /** \copydoc IMemoryOperations::memset */
     void* memset (void* s, int c, size_t n)                 { return _ope.memset (s, c, n);      }

     /** \copydoc IMemoryOperations::memcpy */
     void* memcpy (void* dest, const void* src, size_t n)    { return _ope.memcpy (dest, src, n); }

     /** \copydoc IMemoryOperations::memcmp */
     int   memcmp (const void* s1, const void* s2, size_t n) { return _ope.memcmp (s1, s2, n);    }

     /** \copydoc IMemory::getNbBlocks */
     size_t getNbBlocks () { return 0; }

     /** \copydoc IMemory::getCurrentUsage */
     TotalSize_t getCurrentUsage () { return 0; }

     /** \copydoc IMemory::getMaximumUsage */
     TotalSize_t getMaximumUsage () { return 0; }

protected:

     IMemoryAllocator&  _alloc;
     IMemoryOperations& _ope;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_SYSTEM_IMPL_MEMORY_COMMON_HPP_ */
