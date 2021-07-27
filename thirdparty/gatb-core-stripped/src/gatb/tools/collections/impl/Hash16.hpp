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

/** \file Hash16.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Hash function
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_HASH16_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_HASH16_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Container.hpp>
#include <gatb/tools/collections/api/Bag.hpp>
#include <gatb/tools/misc/impl/Pool.hpp>
#include <gatb/system/impl/System.hpp>

#include <set>
#include <algorithm>
#include <cmath> // for log2f

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
namespace impl          {
/********************************************************************************/


inline uint64_t hash1(uint64_t key, uint64_t seed = 0)
{
    key = (~key) + (key << 21);  // key *= (1 << 21) - 1; key -= 1;
    key = key ^ (key >> 24);
    key = key + (key << 3) + (key << 8);  // key *= 1 + (1 << 3) + (1 << 8)
    key = key ^ (key >> 14);
    key = key + (key << 2) + (key << 4);  // key *= 1 + (1 << 2) + (1 << 4)
    key = key ^ (key >> 28);
    key = key + (key << 31);  // key *= 1 + (1 << 31)
    return key;
}

/** \brief Class providing hash table service with a given memory max usage.
 */
template <typename Item, typename value_type=int> class Hash16
{
	
public :
	
	//shortcut
	typedef misc::impl::cell_ptr_t cell_ptr_t;
	
	typedef struct
	{
		Item graine;
		cell_ptr_t  suiv;
		value_type val;
	} cell;
	
protected:


	
    cell_ptr_t * datah;

    misc::impl::Pool<cell> storage;  // was Item,value_type
    u_int64_t mask ;

    u_int64_t tai;
    u_int64_t nb_elem;
    u_int64_t max_nb_elem;

    /** Shortcut */
    system::IMemory& _memory;

public:


	
    /** Constructor.
     * \param[in] sizeMB : approx max memory to be used by the hash table
     */
    Hash16 (size_t sizeMB) : datah(0), mask(0), tai(0), nb_elem(0), max_nb_elem(0), _memory(system::impl::System::memory())
    {
        int tai_Hash16 = std::max (
            (size_t) ceilf (log2f ((0.1*sizeMB*1024L*1024L)/sizeof(cell_ptr_t))),
            (size_t)1
        );

        /** We check that the provided size is ok. */
        if (tai_Hash16 > 32)  {  throw system::Exception ("Hash16: max size for this hash is 2^32, but ask for %d", tai_Hash16);   }

        tai         = (1LL << tai_Hash16);
        mask        = tai-1 ;
        max_nb_elem = (u_int64_t) (0.8*sizeMB*1024LL*1024LL /sizeof(cell)); // indicative , irrelevant
        //datah       = (cell_ptr_t *) _memory.malloc( tai * sizeof(cell_ptr_t));
		//GR: bug for large values because malloc takes  BlockSize_t (u_int32_t)
		//switching to calloc to avoid problem temporarily, but BlockSize_t should be changed to u_int64_t ?
 		datah       = (cell_ptr_t *) _memory.calloc( tai , sizeof(cell_ptr_t));  //create hashtable

		//printf("Hash16 size asked in MB %zu  tai_Hash16 %i  nb entries %llu \n",sizeMB,tai_Hash16,tai);
		//cell pcell;
		//printf("Hash 16 cell %lli   graine %i suiv %i val %i\n",sizeof(cell),sizeof(pcell.graine),sizeof(pcell.suiv),sizeof(pcell.val));

        _memory.memset (datah,0, tai * sizeof(cell_ptr_t));
    }

	/** Constructor with directly number of entries wished, return really created in nb_created
	 * \param[in] nb_entries : number of entries.
     * \param[in] nb_created : number of created items.
	 */
	Hash16 (u_int64_t nb_entries, u_int64_t * nb_created) : datah(0), mask(0), tai(0), nb_elem(0), max_nb_elem(0), _memory(system::impl::System::memory())
    {
        int tai_Hash16 = std::max (
								   (size_t) ceilf (log2f (nb_entries)),
								   (size_t)1
								   );
		
        /** We check that the provided size is ok. */
        if (tai_Hash16 > 32)  {  throw system::Exception ("Hash16: max size for this hash is 2^32, but ask for %d", tai_Hash16);   }
		
        tai         = (1LL << tai_Hash16);
        mask        = tai-1 ;
        max_nb_elem = (u_int64_t) (10 * tai); // indicative

		if(nb_created!= NULL)
			*nb_created = tai;
 		datah       = (cell_ptr_t *) _memory.calloc( tai , sizeof(cell_ptr_t));  //create hashtable
		
		//printf("Hash16 size asked in MB %zu  tai_Hash16 %i  nb entries %llu \n",sizeMB,tai_Hash16,tai);
		
        _memory.memset (datah,0, tai * sizeof(cell_ptr_t));
    }
	
	u_int64_t getByteSize()
	{
		return storage.getByteSize();
	}
	
    /** Destructor */
    ~Hash16()
    {
        _memory.free(datah);
    }

    /** Clear the content of the hash table. */
    void clear ()
    {
        storage.clear ();
        nb_elem=0;
        _memory.memset (datah,0, tai * sizeof(cell_ptr_t));
    }

    /** Insert an item into the hash table
     * \param[in] graine : key
     * \param[in] value : value
     */
    void insert (Item graine, value_type value)
    {
        unsigned int clef ;
        cell* cell_ptr    = 0;
        cell* newcell_ptr = 0;

        cell_ptr_t  newcell_internal_ptr;

        clef = (unsigned int) (hash1(graine,0) & mask);

        cell_ptr = storage.internal_ptr_to_cell_pointer (datah[clef]);

        while(cell_ptr != NULL &&  cell_ptr->graine != graine)
        {
            cell_ptr = storage.internal_ptr_to_cell_pointer (cell_ptr->suiv);
        }
        if (cell_ptr==NULL) //graine non trouvee , insertion au debut
        {
            newcell_internal_ptr = storage.allocate_cell();
            newcell_ptr          = storage.internal_ptr_to_cell_pointer(newcell_internal_ptr);
            newcell_ptr->val     = value;
            newcell_ptr->graine  = graine;
            newcell_ptr->suiv    = datah[clef];
            datah[clef]          = newcell_internal_ptr;
            nb_elem++;
        }
        else
        {
            cell_ptr->val=value;  // graine trouvee
        }
    }

    /** Insert an item into the hash table
     * \param[in] graine : key
     */
    void insert (Item graine)
    {
        unsigned int clef ;
        cell* cell_ptr, *newcell_ptr;
        cell_ptr_t  newcell_internal_ptr;

        clef = (unsigned int) hash1(graine,0) & mask;

        cell_ptr = storage.internal_ptr_to_cell_pointer(datah[clef]);

        while(cell_ptr != NULL &&  cell_ptr->graine != graine)
        {
            cell_ptr = storage.internal_ptr_to_cell_pointer (cell_ptr->suiv);
        }

        if (cell_ptr==NULL) //graine non trouvee , insertion au debut
        {
            newcell_internal_ptr = storage.allocate_cell();

            newcell_ptr         = storage.internal_ptr_to_cell_pointer(newcell_internal_ptr);
            newcell_ptr->val    = 1;
            newcell_ptr->graine = graine;
            newcell_ptr->suiv   = datah[clef];

            datah[clef] = newcell_internal_ptr;

            nb_elem++;
        }
        else
        {
            (cell_ptr->val)++;  // graine trouvee
        }
    }

	
	static bool sortByKey(const cell &lhs, const cell &rhs) { return lhs.graine < rhs.graine; }
	
	/** Get an iterator for the hash table.
	 * \param[in] sorted : if true, items are iterated in a sorted way (warning: reorder in place so cant acces hash after that !)
	 * \return an iterator over the items of the hash table.
	 */
	//dp::Iterator < std::pair<Item,value_type> >* iterator (bool sorted=false)
	//just get the underlying pool iterator which is simply iteration over multiple arrays, no need to traverse linked list
	dp::Iterator < cell >* iterator (bool sorted=false)
	{
		if(sorted)
		{
			return storage.iteratorsorted(sortByKey);
		}
		else
		{
			return storage.iterator();
		}
	}

	
    /** Get the value for a given key
     * \param[in] graine : key
     * \param[out] val : value to be retrieved
     * \return 1 if the key exists, 0 otherwise.
     */
    int get (Item graine, value_type * val)
    {
        unsigned int clef ;
        cell* cell_ptr;

        clef = (unsigned int) hash1 (graine,0) & mask;

        cell_ptr = storage.internal_ptr_to_cell_pointer(datah[clef]);
        while(cell_ptr != NULL &&  cell_ptr->graine != graine)
        {
            cell_ptr = storage.internal_ptr_to_cell_pointer(cell_ptr->suiv);
        }

        if (cell_ptr==NULL)
        {
            return 0;
        }
        else
        {
            if (val != NULL)  {  *val = cell_ptr->val;  }
            return 1;
        }
    }

    /** Tells whether or not the given key exists in the hash table.
     * \param[in] graine : key
     * \return true if the key exists, 0 otherwise.
     */
    bool contains (Item graine)
    {
        return get (graine,NULL);
    }

    /** Remove a key from the hash table.
     * \param[in] graine : key to be removed from the hash table
     * \param[out] val : value of the key
     * \return true if the key was found, 0 otherwise
     */
    int remove (Item graine, value_type * val)
    {
        unsigned int clef ;
        cell* cell_ptr;
        cell_ptr_t * cellprec_ptr;

        clef = (unsigned int) hash1 (graine,0) & mask;

        cell_ptr = storage.internal_ptr_to_cell_pointer(datah[clef]);
        cellprec_ptr = & (datah[clef]);

        while(cell_ptr != NULL &&  cell_ptr->graine != graine)
        {
            cellprec_ptr = & (cell_ptr->suiv);
            cell_ptr = storage.internal_ptr_to_cell_pointer(cell_ptr->suiv);
        }

        if (cell_ptr==NULL)
        {
            if (val != NULL)  {   *val = 0;  }
            return 0;

        }
        else
        {
            if (val != NULL)  {    *val = cell_ptr->val;  }
            //delete the cell :
            *cellprec_ptr = cell_ptr->suiv ;

            return 1;
        }
    }

    /** Get the number of items in the hash table
     * \return the number of items. */
    u_int64_t size ()  { return nb_elem; }

    /** Get the max number of items allowed by the hash table
     * \return the max number of items. */
    u_int64_t getMaxNbItems ()  { return max_nb_elem; }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_HASH16_HPP_ */
