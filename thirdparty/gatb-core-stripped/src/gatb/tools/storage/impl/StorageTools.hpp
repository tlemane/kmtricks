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

/** \file StorageTools.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Collection interface
 *
 *  This file holds tools about the storage feature
 */

#ifndef _GATB_CORE_TOOLS_STORAGE_IMPL_STORAGE_TOOLS_HPP_
#define _GATB_CORE_TOOLS_STORAGE_IMPL_STORAGE_TOOLS_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Collection.hpp>
#include <gatb/tools/math/NativeInt8.hpp>

#include <gatb/tools/storage/impl/Storage.hpp>
#include <gatb/tools/collections/impl/Bloom.hpp>
#include <gatb/tools/collections/impl/ContainerSet.hpp>


/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace storage   {
namespace impl      {
/********************************************************************************/

/** \brief  Helper class for the storage feature.
 */
class StorageTools
{
public:

    /* why is there code in headers? because else we'd have to instantiate the template, see http://stackoverflow.com/questions/8752837/undefined-reference-to-template-class-constructor */

    /** Singleton method.
     * \return the singleton.
     */
    static StorageTools& singleton() { static StorageTools instance; return instance; }

    /** Save a Collection instance into a group
     * \param[in] group : group where the collection has to be saved
     * \param[in] name : name of the collection saved in the group
     * \param[in] collection : Collection instance to be saved.
     */
    template<typename T>  void saveContainer (Group& group, const std::string& name, collections::Collection<T>* collection)
    {
        collections::Collection<T>* storageCollection = & group.getCollection<T> (name);

        tools::dp::Iterator<T>* it = collection->iterator();   LOCAL(it);
        for (it->first(); !it->isDone(); it->next())  {  storageCollection->insert (it->item());  }
        storageCollection->flush ();
    }

    /** Load a Collection instance from a group
     * \param[in] group : group where the collection has to be load
     * \param[in] name : name of the collection the group
     * \return a Collection instance, loaded from the group
     */
    template<typename T>  collections::Container<T>*  loadContainer (Group& group, const std::string& name)
    {
        collections::Collection<T>*  storageCollection = & group.getCollection<T> (name);
        return new collections::impl::ContainerSet<T> (storageCollection->iterator());
    }

    /** Save a Bloom filter into a group
     * \param[in] group : group where the IBloom instance has to be saved
     * \param[in] name : name of the Bloom filter in the group
     * \param[in] bloom : Bloom filter to be saved
     * \param[in] kmerSize : kmer size (uggly but needed...)
     */
    template<typename T>  void saveBloom (Group& group, const std::string& name, collections::impl::IBloom<T>* bloom, size_t kmerSize)
    {
        collections::Collection<math::NativeInt8>* bloomCollection = & group.getCollection<math::NativeInt8> (name);

        if (bloomMode == 0)
        {
            bloomCollection->insert ((math::NativeInt8*)bloom->getArray(), bloom->getSize());
        }
        else
        {
            // Now we declare an input stream on the collection
            tools::storage::impl::Storage::ostream os (group, name);

            // We write some information in this stream
            os.write (reinterpret_cast<char const*>(bloom->getArray()), bloom->getSize()*sizeof(char));

            // We have to flush the stream in order to be sure everything is ok
            os.flush();
        }

        std::stringstream ss1;  ss1 <<  bloom->getBitSize();
        std::stringstream ss2;  ss2 <<  bloom->getNbHash();
        std::stringstream ss3;  ss3 <<  kmerSize;

        bloomCollection->addProperty ("size",      ss1.str());
        bloomCollection->addProperty ("nb_hash",   ss2.str());
        bloomCollection->addProperty ("type",      bloom->getName());
        bloomCollection->addProperty ("kmer_size", ss3.str());
        bloomCollection->flush (); // R: wasn't there before but I guess this can't hurt
    }

    /** Load a Bloom filter from a group
     * \param[in] group : group where the Bloom filter is
     * \param[in] name : name of the Bloom filter in the group
     * \return the Bloom filter as an instance of IBloom
     */
    template<typename T>  collections::impl::IBloom<T>*  loadBloom (Group& group, const std::string& name)
    {
        /** We retrieve the raw data buffer for the Bloom filter. */
        tools::collections::Collection<tools::math::NativeInt8>* bloomArray = & group.getCollection<tools::math::NativeInt8> (name);

        /** We create the Bloom fiter. */
        tools::collections::impl::IBloom<T>* bloom = tools::collections::impl::BloomFactory::singleton().createBloom<T> (
            bloomArray->getProperty("type"),
            bloomArray->getProperty("size"),
            bloomArray->getProperty("nb_hash"),
            bloomArray->getProperty("kmer_size")
        );

        if (bloomMode == 0)
        {
            /** We set the bloom with the provided array given as an iterable of NativeInt8 objects. */
            bloomArray->getItems ((tools::math::NativeInt8*&)bloom->getArray());
        }
        else
        {
            tools::storage::impl::Storage::istream is (group, name);

            // We read the data from the input stream
            is.read (reinterpret_cast<char*>(bloom->getArray()), bloom->getSize()*sizeof(char));
        }

        /** We return the result. */
        return bloom;
    }

private:

    /** We keep the possibility to load/save Bloom filters in two different ways.
     * The old one (with 'insert') has the drawback that the Bloom filter was read/written in
     * one shot, which made big memory usage by h d f 5 (one buffer for memory, one buffer for file).
     * The new way uses io streams for h d f 5 (splits the read/write operations in chunks), which
     * may use less memory. */
    static const int bloomMode = 1;  // 0: uses insert    1: uses Storage:[oi]stream
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_STORAGE_IMPL_STORAGE_TOOLS_HPP_ */
