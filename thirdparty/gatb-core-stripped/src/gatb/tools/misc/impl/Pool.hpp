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

/** \file Pool.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Tool for pooling objects
 */

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_POOL_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_POOL_HPP_

/********************************************************************************/

#include <gatb/system/impl/System.hpp>
#include <queue>          // std::priority_queue

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/* Cette class dÈfinit une pool memoire pour allocation rapide de la table de hachage
 * utilisee quand seed >14 */
	
//now pass as template the cell type
	
//had to move this outside the class (otherwise recursive def of cell template that may contain cell_ptr_t ... )
typedef u_int32_t  cell_ptr_t;

	
template <typename cell>  class Pool
//template <typename graine_type, typename value_type=int>  class Pool
{
public:

	/** Default constructor.
     * \param[in] tai :  2^20  1 M cells *16 o    blocs de 16 Mo
     * \param[in] N : 2^12  soit 4 G cells max
     * */
    Pool (size_t tai=1048576, size_t N=4096) : TAI_POOL(tai), N_POOL(N)
    {
        n_pools = 0; n_cells=0;
        //allocation table de pool :
        tab_pool = (cell**)  MALLOC (N_POOL*sizeof(cell*) );

        tab_pool[0]=0; n_pools++; // la premiere pool est NULL, pour conversion null_internal -> null

        //allocation de la premiere pool :
        pool_courante =(cell*)  MALLOC (TAI_POOL*sizeof(cell) );
        tab_pool[n_pools] = pool_courante;
        n_pools++;
    }

    /**  Destructeur  */
    ~Pool()
    {
        // la pool 0 est NULL
        for(size_t i=1;i<n_pools;i++)  {  FREE ( tab_pool[i] );  }

        FREE (tab_pool);
    }

	
	u_int64_t getByteSize()
	{
		return ((n_pools-2) * TAI_POOL*sizeof(cell)); // I do not want to count the firs tdefault allocated pool
	}
	
	
    /**  allocate cell, return internal pointer type ( 32bits) */
    cell_ptr_t  allocate_cell()
    {
        cell_ptr_t internal_adress = 0;
        // ncells = nb de cells deja utilisees
        if (n_cells <TAI_POOL)
        {
            internal_adress  = n_pools -1;    // low 12 bits  : pool number
            internal_adress |= n_cells << 12; // 20 high bits : cell number in pool
            n_cells ++;

            return internal_adress;
        }
        else // la piscine est pleine, on en alloue une nouvelle
        {
            if (n_pools>= N_POOL)
            {
                // will happen when  4G cells are allocated, representing 64 Go
                throw system::Exception ("Internal memory allocator is full!");
            }
            pool_courante =(cell*)  MALLOC (TAI_POOL*sizeof(cell) );
            tab_pool[n_pools] = pool_courante;
            n_pools++;
            n_cells = 1;

            internal_adress = n_pools -1;
            // 20 high bits are 0

            return internal_adress;
        }
    }

    /** */
    cell*  internal_ptr_to_cell_pointer(cell_ptr_t internal_ptr)
    {
        unsigned int numpool =  internal_ptr & 4095;
        unsigned int numcell =  internal_ptr >> 12;

        return (tab_pool[numpool] + numcell);
    }

    /** vide toutes piscines (garde juste une pool vide)  */
    void  clear ()
    {
        for(size_t i=2;i<n_pools;i++) // garde la premiere pool pour usage futur
        {
            FREE ( tab_pool[i] );
        }
		memset(tab_pool[1],0,TAI_POOL*sizeof(cell));
		
        //on repasse sur premiere pool
        pool_courante = tab_pool[1];
        n_cells=0;
        n_pools=2;
    }


	
	//sort the pools according to some comparator
	//warning this will reorder cells and thus  making existing pointers to cells irrelevant
	//but useful for  e.g. sirted iterator of cells
	template <typename Comparator>
	void sortPools(Comparator comparator)
	{
		// les pool pleines
		for(size_t i=1;i<(n_pools-1);i++)
		{
			std::sort( tab_pool[i],  tab_pool[i]  + TAI_POOL, comparator);
		}
		
		// la pool en cours de remplissage
		std::sort( tab_pool[n_pools-1],  tab_pool[n_pools-1]  + n_cells, comparator);

	}
	
	
	////////simple iterator over all cells
	template <typename Comparator>
	dp::Iterator < cell >* iteratorsorted (Comparator comparator)
	{
		//first sort each pool with std sort
		 this->sortPools(comparator);

		//then iterate with a merge sort
		return new IteratorSorted(*this);
	}
	
	
	//todo template also this with a comparator
	class IteratorSorted : public tools::dp::Iterator < cell >
	{
	
	public:
		typedef std::pair<int, cell *> cellpair_t; //id pointer of pool , cell *

		struct sortcellpair { bool operator() (cellpair_t &l,cellpair_t &r) { return !(  (* l.second).graine <=  (* r.second).graine );  }  } ;

		IteratorSorted (Pool<cell>& aRef) : ref(aRef), done(true)  {}
		
		/** \copydoc tools::dp::Iterator::first */
		void first()
		{

			for(size_t i=1;i< ref.n_pools;i++)
			{
				pq.push( cellpair_t(i,  (cell *)  &(ref.tab_pool[i][0])  )   );
			}
			next();
		}
		
		/** \copydoc tools::dp::Iterator::next */
		void next()
		{
			
			done =  (pq.size() == 0);
			
			if(!done)
			{
				cellpair_t current_pair = pq.top() ; pq.pop();
				*this->_item = * (current_pair.second);
				

				
				//push the next cell of this list  if any
				unsigned int cell_number =  current_pair.second   -  ref.tab_pool[current_pair.first] ;
				unsigned int current_pool = current_pair.first;
				if( (current_pool < (ref.n_pools -1)) &&  ((cell_number+1)  < ref.TAI_POOL)  ) // inside a full pool, and cells remaining
				{
					pq.push( cellpair_t(current_pool, & (ref.tab_pool[current_pool][cell_number+1])  )   );
				}
				else if ( (current_pool == (ref.n_pools -1)) && ((cell_number+1) < ref.n_cells) ) // inside last pool, and cells remaining
				{
					pq.push( cellpair_t(current_pool, & (ref.tab_pool[current_pool][cell_number+1])  )   );
				}
				//otherwise at end of array, dont push anything
			}
		}
		
		/** \copydoc tools::dp::Iterator::isDone */
		bool isDone ()   {  return done; }
		
		/** \copydoc tools::dp::Iterator::item */
		cell& item ()     { return *this->_item; }
		
	private:
		std::priority_queue< cellpair_t, std::vector<cellpair_t>,   sortcellpair > pq;
		Pool<cell>&  ref;
		bool         done;
		
	};
	
	//////
	
	
	
	
	
	
	////////simple iterator over all cells
	dp::Iterator < cell >* iterator ()
	{
		return new Iterator(*this);
	}
	
	
	/************************************************************/
	//avec std::pair ? pour avoir Item, value_type
	class Iterator : public tools::dp::Iterator <  cell  >
	{
	public:
		
		Iterator (Pool<cell>& aRef) : ref(aRef), done(true)  {}
		
		/** \copydoc tools::dp::Iterator::first */
		void first()
		{
			_current_pool = 1; // first pool
			_current_cell = 0;
			done        = ref.n_cells < 1;
			if(!done)
				*this->_item = ref.tab_pool[_current_pool][_current_cell];
			
			_current_cell++; // next cell that should be read
		}
		
		/** \copydoc tools::dp::Iterator::next */
		void next()
		{

			if(_current_pool < (ref.n_pools -1) && _current_cell  < ref.TAI_POOL ) // inside a full pool, and cells remaining
			{
				*this->_item = ref.tab_pool[_current_pool][_current_cell];
				_current_cell++;
				return;
			}
			else if (_current_pool < (ref.n_pools -1) && _current_cell == ref.TAI_POOL ) // inside a full pool but no cells remaining
			{
				//go to next pool
				_current_pool++;
				_current_cell = 0;
				*this->_item = ref.tab_pool[_current_pool][_current_cell];
				_current_cell++;
			}
			else if (_current_pool == (ref.n_pools -1) && _current_cell < ref.n_cells) // in last pool and cells remaining
			{
				*this->_item = ref.tab_pool[_current_pool][_current_cell];
				_current_cell++;
			}
			else // last pools and no cells left, done
			{
				done = true;
			}
		}
		
		/** \copydoc tools::dp::Iterator::isDone */
		bool isDone ()   {  return done; }
		
		/** \copydoc tools::dp::Iterator::item */
		cell& item ()     { return *this->_item; }
		
	private:
		
		unsigned int _current_pool;
		unsigned int _current_cell;

		Pool<cell>&  ref;
		
		bool         done;
	};
	

	
private:

    /** table de cell, pour usage courant */
    cell* pool_courante;

    /** stockage de tous les pointeurs pool */
    cell** tab_pool;

    /** nombre de piscines remplies */
    unsigned int n_pools;

    /**  niveau de remplissage de la piscine courante */
    unsigned int n_cells;

    size_t TAI_POOL;
    size_t N_POOL;
};

/********************************************************************************/

//make it an allocator usable by std vector ?
class MemAllocator
{
public:

    //clear all previous allocs, and alloc pool capacity
    void reserve(u_int64_t size)
    {
        if(mainbuffer !=NULL)
        {
            FREE (mainbuffer);
        }

        /** We add a little bit of memory in case "align" method is called often.
         * We allow max alignment of 16 bytes per core, plus some extra memory. */
        size_t extraMem = 16*_nbCores + 1024;

        capacity   = size+extraMem;
        mainbuffer = (char*) CALLOC(capacity,1);
        used_space = 0;
    }

    //should be thread safe
    char* pool_malloc(u_int64_t requested_size, const char* message="")
    {
        u_int64_t synced_used_space = __sync_fetch_and_add(&used_space, requested_size);

        if (requested_size > (capacity - synced_used_space))
        {
            __sync_fetch_and_add(&used_space, -requested_size);

            std::cerr << "Allocation failed for " << std::to_string(requested_size) << ". Current is " << std::to_string(used_space) << ", capacity is " << std::to_string(capacity);
            std::cerr << ". Partition -> " << message << std::endl;
            throw system::Exception ("Pool allocation failed for %lld bytes (%s). Current usage is %lld and capacity is %lld",
                requested_size, message, used_space, capacity
            );

            return NULL;
        }

        if (mainbuffer == NULL)
        {

            std::cerr << "Allocation failed, main buffer is null. Current is " << std::to_string(used_space) << ", capacity is " << std::to_string(capacity) << std::endl;
            throw system::Exception ("Pool allocation failed for %lld bytes (%s), mainbuffer is null?. Current usage is %lld and capacity is %lld",
                requested_size, message, used_space, capacity
                );
        }
        return mainbuffer + synced_used_space;
    }

    /** Force alignment. */
    void align (u_int8_t alignBytes)
    {
        size_t offset = alignBytes-1 + sizeof(char*);
        char* current = mainbuffer + used_space;
        char* buffer  =  (char*)(((size_t)(current)+offset)&~(alignBytes-1));

        /** We update the used space. */
        used_space = (u_int64_t)(buffer-mainbuffer);
    }

    u_int64_t getCapacity ()  {  return capacity;   }

    u_int64_t getUsedSpace()  {  return used_space; }


    void free_all()
    {
        used_space = 0;
    }

    MemAllocator(size_t nbCores=0) : mainbuffer(NULL),capacity(0),used_space(0), _nbCores(nbCores), _synchro(0)
    {
        setSynchro (system::impl::System::thread().newSynchronizer());
    }

    ~MemAllocator()
    {
        if (mainbuffer != NULL)  {  FREE (mainbuffer);  }
        setSynchro (0);
    }

    system::ISynchronizer* getSynchro()  { return _synchro; }

private :
    char*     mainbuffer;
    u_int64_t capacity; //in bytes
    u_int64_t used_space;

    size_t _nbCores;

    system::ISynchronizer* _synchro;
    void setSynchro (system::ISynchronizer* synchro) { SP_SETATTR(synchro); }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_POOL_HPP_ */
