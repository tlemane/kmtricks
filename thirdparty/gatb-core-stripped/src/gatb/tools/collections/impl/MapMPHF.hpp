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

/** \file MapMPHF.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Minimal Perfect Hash Function
 */

#ifndef _GATB_CORE_TOOLS_COLLECTION_MAP_MPHF_HPP_
#define _GATB_CORE_TOOLS_COLLECTION_MAP_MPHF_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/collections/impl/BooPHF.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <vector>

/********************************************************************************/
namespace gatb        {
	namespace core        {
		namespace tools       {
			namespace collections {
				namespace impl        {
					/********************************************************************************/
					
					/** \brief hash table implementation
					 *
					 * This hash table implementation uses a minimal perfect hash function (MPHF) for
					 * identifying the keys with a unique number in [0..N-1] where N is the number of items.
					 *
					 * Using BooPHF, the memory usage is about 3-4 bits per key.
					 *
					 * The values can be stored in a simple vector. The keys are not stored in memory, only
					 * the mphf is needed.
					 *
					 * Note that such an implementation can't afford to add items into the map (it's static).
					 */
					template <class Key, class Value, class Adaptator=AdaptatorDefault<Key> >
					class MapMPHF : public system::SmartPointer
					{
					public:
						
						/** Hash type. */
						typedef BooPHF<Key, Adaptator> Hash;
						
						/** Default constructor. */
						MapMPHF () : hash() {}
						
						/** Build the hash function from a set of items.
						 * \param[in] keys : iterable over the keys of the hash table
						 * \param[in] progress : listener called during the building of the MPHF
						 */
						void build (tools::collections::Iterable<Key>& keys, int nbThreads = 1, tools::dp::IteratorListener* progress=0)
						{
							/** We build the hash function. */
							hash.build (&keys, nbThreads, progress);
							
							/** We resize the vector of Value objects. */
							data.resize (keys.getNbItems());
							clearData();
							initDiscretizationScheme();
						}
						
						
						
						// discretization scheme to store abundance values from 0 to 50000 on 8 bits
						// with  5% error maximum
						// from 0     to 70     :   step = 1              (70 buckets)
						// from 70    to 100    :   step = 2              (15 buckets)
						// from 100   to 500    :   step = 10             (40 buckets)
						// from 500   to 1000   :   step = 20             (25 buckets)
						// from 1000  to 5000   :   step = 100            (40 buckets)
						// from 5000  to 10000  :   step = 200            (25 buckets)
						// from 10000 to 50000  :   step = 1000           (40 buckets)
						//
						//to change discretization scheme, change the values in _abundanceDiscretization below
					
						void initDiscretizationScheme()
						{
							_abundanceDiscretization.resize(257);
							
							int total =0;
							_abundanceDiscretization[0] = 0;
							int idx=1;
							for(int ii=1; ii<= 70; ii++,idx++ )
							{
								total += 1;
								_abundanceDiscretization[idx] = total ;
							}
							
							for(int ii=1; ii<= 15; ii++,idx++ )
							{
								total += 2;
								_abundanceDiscretization[idx] = total  ;
							}
							
							for(int ii=1; ii<= 40; ii++,idx++ )
							{
								total += 10;
								_abundanceDiscretization[idx] = total  ;
							}
							for(int ii=1; ii<= 25; ii++,idx++ )
							{
								total += 20;
								_abundanceDiscretization[idx] = total  ;
							}
							for(int ii=1; ii<= 40; ii++,idx++ )
							{
								total += 100;
								_abundanceDiscretization[idx] = total  ;
							}
							for(int ii=1; ii<= 25; ii++,idx++ )
							{
								total += 200;
								_abundanceDiscretization[idx] = total  ;
							}
							for(int ii=1; ii<= 40; ii++,idx++ )
							{
								total += 1000;
								_abundanceDiscretization[idx] = total  ;
							}
							//
							
							
							
							_abundanceDiscretization[256] = total;

							
							
							
							
							/*
							for(int ii=0; ii<_abundanceDiscretization.size(); ii++ )
							{
								printf("disc[%i] = %i \n", ii,_abundanceDiscretization[ii] );
							}
							*/
							
						}
						
						/* use the hash from another MapMPHF class. hmm is this smartpointer legit?
						 * also allocate n/x data elements
						 */
						void useHashFrom (MapMPHF *other, int x = 1)
						{
							hash = other->hash;
							
							/** We resize the vector of Value objects. */
							data.resize ((unsigned long)((hash.size()) / (unsigned long)x) + 1LL); // that +1 and not (hash.size+x-1) / x
							
							clearData();
						}
						
						/** Save the hash function into a Group object.
						 * \param[out] group : group where to save the MPHF
						 * \param[in] name : name of the saved MPHF
						 * \return the number of bytes of the saved data.
						 */
						size_t save (tools::storage::impl::Group& group, const std::string& name)  {  return hash.save (group, name);  }
						
						/** Load hash function from a Group
						 * \param[in] group : group where to load the MPHF from
						 * \param[in] name : name of the MPHF
						 */
						void load (tools::storage::impl::Group& group, const std::string& name)
						{
							/** We load the hash function. */
							size_t nbKeys = hash.load (group, name);
							
							/** We resize the vector of Value objects. */
							data.resize (nbKeys);
							clearData();
							initDiscretizationScheme();
						}
						
						/** Get the value for a given key
						 * \param[in] key : the key
						 * \return the value associated to the key. */
						Value& operator[] (const Key& key)  {
							return data[hash(key)];
						}
						
						/** Get the value for a given index
						 * \param[in] code : the key
						 * \return the value associated to the key. */
						Value& at (typename Hash::Code code)  {
							return data[code];
						
						}
						
						Value& at (const Key& key)  {
							return data[hash(key)];
						}
						
                        int abundanceAt (const Key& key)  {
							return floorf((_abundanceDiscretization [data[hash(key)]]  +  _abundanceDiscretization [data[hash(key)]+1])/2.0);
						}
	
                        int abundanceAt (typename Hash::Code code)  {
							return floorf((_abundanceDiscretization [data[code]]  +  _abundanceDiscretization [data[code]+1])/2.0);
						}
						
						/** Get the hash code of the given key. */
						typename Hash::Code getCode (const Key& key) { return hash(key); }
						
						/** Get the number of keys.
						 * \return keys number. */
						size_t size() const { return hash.size(); }
						
						void clearData() { 
							for (unsigned long i = 0; i < data.size(); i ++)
								data[i] = 0;
						}
						
						std::vector<int>   _abundanceDiscretization;

					private:
						
						Hash               hash;
						std::vector<Value> data;
						
						
					};
					
					/********************************************************************************/
				} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTION_MAP_MPHF_HPP_ */
