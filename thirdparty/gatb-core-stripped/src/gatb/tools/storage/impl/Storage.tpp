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

/********************************************************************************/
namespace gatb  {  namespace core  {  namespace tools  {  namespace storage  {  namespace impl {
/********************************************************************************/

/********************************************************************************
     #####   #######  #        #        #######   #####   #######
    #     #  #     #  #        #        #        #     #     #
    #        #     #  #        #        #        #           #
    #        #     #  #        #        #####    #           #
    #        #     #  #        #        #        #           #
    #     #  #     #  #        #        #        #     #     #      ##
     #####   #######  #######  #######  #######   #####      #      ##
********************************************************************************/

/*********************************************************************
*********************************************************************/
template <class Item>
inline CollectionNode<Item>::CollectionNode (
    StorageFactory* factory, 
    ICell* parent, 
    const std::string& id, 
    collections::Collection<Item>* ref
)
    : Cell(parent,id), collections::impl::CollectionAbstract<Item> (ref->bag(), ref->iterable()), _factory(factory), _ref(0)
{
    /** We get a token on the referred Collection instance. */
    setRef(ref);
}

/*********************************************************************
*********************************************************************/
template <class Item>
inline CollectionNode<Item>::~CollectionNode()
{
    /** We release the token of the Collection instance. */
    setRef(0);
}

/*********************************************************************
*********************************************************************/
template <class Item>
inline void CollectionNode<Item>::remove ()
{
    /** We delegate the job to the referred Collection instance. */
    _ref->remove();
}

/*********************************************************************
*********************************************************************/
template <class Item>
inline void CollectionNode<Item>::addProperty (const std::string& key, const std::string value)  
{  
    _ref->addProperty (key, value); 
}

/*********************************************************************
*********************************************************************/
template <class Item>
inline std::string  CollectionNode<Item>::getProperty (const std::string& key) 
{  
    return _ref->getProperty (key); 
}

/*********************************************************************
*********************************************************************/
template <class Item>
inline collections::Collection<Item>* CollectionNode<Item>::getRef ()  
{ 
    return _ref; 
}

/**********************************************************************
             #####   ######   #######  #     #  ######
            #     #  #     #  #     #  #     #  #     #
            #        #     #  #     #  #     #  #     #
            #  ####  ######   #     #  #     #  ######
            #     #  #   #    #     #  #     #  #
            #     #  #    #   #     #  #     #  #
             #####   #     #  #######   #####   #
**********************************************************************/

/*********************************************************************
*********************************************************************/
inline Group::Group (StorageFactory* factory, ICell* parent, const std::string& name) 
    : Cell(parent, name), _factory(factory) 
{
}

/*********************************************************************
*********************************************************************/
inline Group::~Group()
{
    /** We release the collections. */
    for (size_t i=0; i<_collections.size(); i++)  {  _collections[i]->forget();  }

    /** We release the partitions. */
    for (size_t i=0; i<_partitions.size(); i++)  { _partitions[i]->forget(); }

    /** We release the groups. */
    for (size_t i=0; i<_groups.size(); i++)  { _groups[i]->forget(); }
}

/*********************************************************************
*********************************************************************/
inline Group& Group::getGroup (const std::string& name)
{
    Group* group=0;  for (size_t i=0; !group && i<_groups.size(); i++)  {  if (_groups[i]->getId() == name)  { group = _groups[i]; }  }

    if (group == 0)
    {
        group = _factory->createGroup (this, name);
        _groups.push_back (group);
        group->use ();
    }
    return *group;
}

/*********************************************************************
*********************************************************************/
template <class Type>  
inline Partition<Type>&  Group::getPartition (const std::string& name, size_t nb)
{
    Partition<Type>* result = _factory->createPartition<Type> (this, name, nb);
    _partitions.push_back(result);
    result->use();
    return *result;
}

/*********************************************************************
*********************************************************************/
template <class Type>  
inline CollectionNode<Type>&  Group::getCollection (const std::string& name)
{
    CollectionNode<Type>* result = _factory->createCollection<Type> (this, name, 0);
    _collections.push_back (result);
    result->use ();
    return *result;
}

/*********************************************************************
*********************************************************************/
inline void Group::remove ()
{
    /** We remove the collections. */
    for (size_t i=0; i<_collections.size(); i++)  {  _collections[i]->remove ();  }

    /** We remove the partitions. */
    for (size_t i=0; i<_partitions.size(); i++)   { _partitions[i]->remove(); }

    /** We remove the children groups. */
    for (size_t i=0; i<_groups.size(); i++)       { _groups[i]->remove(); }
}

/**********************************************************************
######      #     ######   #######  ###  #######  ###  #######  #     #
#     #    # #    #     #     #      #      #      #   #     #  ##    #
#     #   #   #   #     #     #      #      #      #   #     #  # #   #
######   #     #  ######      #      #      #      #   #     #  #  #  #
#        #######  #   #       #      #      #      #   #     #  #   # #
#        #     #  #    #      #      #      #      #   #     #  #    ##
#        #     #  #     #     #     ###     #     ###  #######  #     #
**********************************************************************/

/*********************************************************************
*********************************************************************/
template<typename Type>
inline Partition<Type>::Partition (StorageFactory* factory, ICell* parent, const std::string& id, size_t nbCollections)
    : Group (factory, parent, id), _factory(factory), _typedCollections(nbCollections), _synchro(0)
{
    /** We create a synchronizer to be shared by the collections. */
    _synchro = system::impl::System::thread().newSynchronizer();

    /** We want to instantiate the wanted number of collections. */
    for (size_t i=0; i<_typedCollections.size(); i++)
    {
        /** We define the name of the current partition as a mere number. */
        std::stringstream ss;  ss << i;

        CollectionNode<Type>* result = _factory->createCollection<Type> (this, ss.str(), _synchro);

        /** We add the collection node to the dedicated vector and take a token for it. */
        (_typedCollections [i] = result)->use ();
    }
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline Partition<Type>::~Partition ()
{
    /** We release the token for each collection node. */
    for (size_t i=0; i<_typedCollections.size(); i++)  { _typedCollections[i]->forget ();  }

    delete _synchro;
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline size_t Partition<Type>::size()  const  
{ 
    return _typedCollections.size(); 
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline collections::Collection<Type>& Partition<Type>::operator[] (size_t idx)  
{  
    return  * _typedCollections[idx]->getRef();  
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline dp::Iterator<Type>* Partition<Type>::iterator ()
{
    std::vector <dp::Iterator<Type>*> iterators;
    for (size_t i=0; i<this->size(); i++) { iterators.push_back ((*this)[i].iterator()); }
    return new dp::impl::CompositeIterator<Type> (iterators);
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline int64_t Partition<Type>::getNbItems ()
{
    int64_t result = 0;   for (size_t i=0; i<this->size(); i++) { result += (*this)[i].getNbItems(); }
    //std::cout << "returning getNbItems from " << this->size() << " partitions: " << result << std::endl;
    return result;
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline u_int64_t Partition<Type>::getSizeItems ()
{
    u_int64_t result = 0;   for (size_t i=0; i<this->size(); i++) { result += (*this)[i].getNbItems() * sizeof(Type); }
    return result;
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline int64_t Partition<Type>::estimateNbItems ()
{
    int64_t result = 0;   for (size_t i=0; i<this->size(); i++) { result += (*this)[i].estimateNbItems(); }
    return result;
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline void Partition<Type>::flush ()
{
    for (size_t i=0; i<_typedCollections.size(); i++)  { _typedCollections[i]->flush();  }
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline void Partition<Type>::remove ()
{
    /** We remove each collection of this partition. */
    for (size_t i=0; i<_typedCollections.size(); i++)  { _typedCollections[i]->remove ();  }

    /** We call the remove of the parent class. */
    Group::remove ();
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline PartitionCache<Type>::PartitionCache (Partition<Type>& ref, size_t nbItemsCache, system::ISynchronizer* synchro)
    :  _ref(ref), _nbItemsCache(nbItemsCache), _synchro(synchro), _synchros(ref.size()), _cachedCollections(ref.size())
{
    /** We create the partition files. */
    for (size_t i=0; i<_cachedCollections.size(); i++)
    {
        if(synchro==0)
        {
            _synchros[i] = system::impl::System::thread().newSynchronizer();
        }
        else
        {
            _synchros[i] = synchro;
        }
        _synchros[i]->use();

        _cachedCollections[i] = new collections::impl::CollectionCache<Type> (ref[i], nbItemsCache, _synchros[i]);
        _cachedCollections[i]->use ();
    }
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline PartitionCache<Type>::PartitionCache (const PartitionCache<Type>& p)
    : _ref(p._ref), _nbItemsCache(p._nbItemsCache), _synchro(p._synchro) , _synchros(p._synchros), _cachedCollections(p.size())
{
    /** We create the partition files. */
    for (size_t i=0; i<_cachedCollections.size(); i++)
    {
        PartitionCache<Type>& pp = (PartitionCache<Type>&)p;

        _cachedCollections[i] = new collections::impl::CollectionCache<Type> (pp[i].getRef(), p._nbItemsCache, p._synchros[i]);
        _cachedCollections[i]->use ();
        _synchros[i]->use();

    }
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline PartitionCache<Type>::~PartitionCache ()
{
    flush ();
    for (size_t i=0; i<_cachedCollections.size(); i++)  {
        _cachedCollections[i]->forget ();
        _synchros[i]->forget();
    }
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline size_t  PartitionCache<Type>::size() const 
{ 
    return _cachedCollections.size(); 
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline collections::impl::CollectionCache<Type>&   PartitionCache<Type>::operator[] (size_t idx)  
{  
    return * _cachedCollections[idx];  
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline void PartitionCache<Type>::flush ()   
{  
    for (size_t i=0; i<_cachedCollections.size(); i++)  { _cachedCollections[i]->flush();  }  
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline void PartitionCache<Type>::remove ()  
{  
    for (size_t i=0; i<_cachedCollections.size(); i++)  { _cachedCollections[i]->remove ();  } 
}
    
/*********************************************************************
 *********************************************************************/
template<typename Type>
inline PartitionCacheSorted<Type>::PartitionCacheSorted (Partition<Type>& ref, size_t nbItemsCache, u_int32_t max_memory, system::ISynchronizer* synchro)
:  _ref(ref), _nbItemsCache(nbItemsCache), _synchro(synchro), _cachedCollections(ref.size()) , _synchros(ref.size()) , _outsynchros(ref.size()), _sharedBuffers(ref.size()), _idxShared(ref.size()), _max_memory(max_memory),_numthread(-1),_nbliving(0)
{

    //todo should also take into account the temp buffer used for sorting, ie max nb threads * buffer, so should be
   // _sharedBuffersSize = (_max_memory*system::MBYTE / (_cachedCollections.size()  + nbthreads) ) / sizeof(Type); //in nb elems
    
    _sharedBuffersSize = (_max_memory*system::MBYTE / _cachedCollections.size()  ) / sizeof(Type); //in nb elems
    _sharedBuffersSize = std::max(2*nbItemsCache,_sharedBuffersSize);
   // _sharedBuffersSize = 1*system::MBYTE  / sizeof(Type);

    printf("sort cache of %zu elems \n",_sharedBuffersSize);
    /** We create the partition files. */
    for (size_t i=0; i<_cachedCollections.size(); i++)
    {

        if(synchro==0)
        {
            _synchros[i] = system::impl::System::thread().newSynchronizer();
        }
        else
        {
            _synchros[i] = synchro;
        }
        
        _outsynchros[i] = system::impl::System::thread().newSynchronizer();
        
        _synchros[i]->use();
        _outsynchros[i]->use();
        
       // printf("Creating outsync %p   \n",_outsynchros[i] );

        _sharedBuffers[i] = (Type*) system::impl::System::memory().calloc (_sharedBuffersSize, sizeof(Type));
        //system::impl::System::memory().memset (_sharedBuffers[i], 0, _sharedBuffersSize*sizeof(Type));
       // printf("Creating shared buffer %p   \n",_sharedBuffers[i] );


        
        _cachedCollections[i] = new collections::impl::CollectionCacheSorted<Type> (ref[i], nbItemsCache,_sharedBuffersSize, _synchros[i],_outsynchros[i],_sharedBuffers[i],&_idxShared[i]);
        //
        _cachedCollections[i]->use ();
    }
}

/*********************************************************************
 *********************************************************************/
template<typename Type>
inline PartitionCacheSorted<Type>::PartitionCacheSorted (const PartitionCacheSorted<Type>& p)
: _ref(p._ref), _nbItemsCache(p._nbItemsCache), _synchro(p._synchro), _cachedCollections(p.size()) , _synchros(p._synchros), _outsynchros(p._outsynchros), _sharedBuffers(0),_idxShared(0)
{

   //_numthread =  __sync_fetch_and_add (&p._nbliving, 1);

    /** We create the partition files. */
    for (size_t i=0; i<_cachedCollections.size(); i++)
    {

        _synchros[i]->use();
        _outsynchros[i]->use();
        
        PartitionCacheSorted<Type>& pp = (PartitionCacheSorted<Type>&)p;
        
        _cachedCollections[i] = new collections::impl::CollectionCacheSorted<Type> (pp[i].getRef(), p._nbItemsCache,p._sharedBuffersSize, p._synchros[i],p._outsynchros[i],p._sharedBuffers[i],&(pp._idxShared[i]) ); //
        _cachedCollections[i]->use ();
    }
}
/*********************************************************************
 *********************************************************************/
template<typename Type>
inline PartitionCacheSorted<Type>::~PartitionCacheSorted ()
{
    //printf("destruc parti cache sorted tid %i \n",_numthread);
    flush ();

    for (size_t i=0; i<_cachedCollections.size(); i++)  {
        //destruction of collection will also flush the output bag (which is shared between all threads) , so need to lock it also
        _outsynchros[i]->lock();
        _cachedCollections[i]->forget ();
        _outsynchros[i]->unlock();

        
        _synchros[i]->forget();
        _outsynchros[i]->forget();
    }
    
    for (typename std::vector<Type*>::iterator it = _sharedBuffers.begin() ; it != _sharedBuffers.end(); ++it)
    {
        //printf("destroying shared buffer %p \n",*it);
        system::impl::System::memory().free (*it);
    }

}

/*********************************************************************
 *********************************************************************/
template<typename Type>
inline size_t  PartitionCacheSorted<Type>::size() const
{
    return _cachedCollections.size();
}

/*********************************************************************
 *********************************************************************/
template<typename Type>
inline collections::impl::CollectionCacheSorted<Type>&   PartitionCacheSorted<Type>::operator[] (size_t idx)
{
    return * _cachedCollections[idx];
}

/*********************************************************************
 *********************************************************************/
template<typename Type>
inline void PartitionCacheSorted<Type>::flush ()
{
    //printf("PartitionCacheSorted flush  this %p \n",this);
    for (size_t i=0; i<_cachedCollections.size(); i++)  { _cachedCollections[i]->flush();  }
}

/*********************************************************************
 *********************************************************************/
template<typename Type>
inline void PartitionCacheSorted<Type>::remove ()
{  
    for (size_t i=0; i<_cachedCollections.size(); i++)  { _cachedCollections[i]->remove ();  } 
}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#include <gatb/tools/storage/impl/StorageFile.hpp>

/********************************************************************************/
namespace gatb  {  namespace core  {  namespace tools  {  namespace storage  {  namespace impl {
/********************************************************************************/

/********************************************************************************
        #######     #      #####   #######  #######  ######   #     #
        #          # #    #     #     #     #     #  #     #   #   #
        #         #   #   #           #     #     #  #     #    # #
        #####    #     #  #           #     #     #  ######      #
        #        #######  #           #     #     #  #   #       #
        #        #     #  #     #     #     #     #  #    #      #
        #        #     #   #####      #     #######  #     #     #
********************************************************************************/

/*********************************************************************
*********************************************************************/
inline Storage* StorageFactory::create (const std::string& name, bool deleteIfExist, bool autoRemove, bool dont_add_extension, bool append)
{
    switch (_mode)
    {
        case STORAGE_FILE:  return StorageFileFactory::createStorage (name, deleteIfExist, autoRemove);
        case STORAGE_GZFILE:  return StorageGzFileFactory::createStorage (name, deleteIfExist, autoRemove);
        case STORAGE_COMPRESSED_FILE:  return StorageSortedFactory::createStorage (name, deleteIfExist, autoRemove);
        default:            throw system::Exception ("Unknown mode in StorageFactory::createStorage");
    }
}

/*********************************************************************
*********************************************************************/
inline bool StorageFactory::exists (const std::string& name)
{
    switch (_mode)
    {
        case STORAGE_FILE:              return StorageFileFactory::exists (name);
        case STORAGE_GZFILE:            return StorageGzFileFactory::exists (name);
        case STORAGE_COMPRESSED_FILE:   return StorageSortedFactory::exists (name);
        default:            throw system::Exception ("Unknown mode in StorageFactory::exists");
    }
}

/*********************************************************************
*********************************************************************/
inline Group* StorageFactory::createGroup (ICell* parent, const std::string& name)
{
    switch (_mode)
    {
        case STORAGE_FILE:  return StorageFileFactory::createGroup (parent, name);
        case STORAGE_GZFILE:  return StorageGzFileFactory::createGroup (parent, name);
        case STORAGE_COMPRESSED_FILE:  return StorageSortedFactory::createGroup (parent, name);

        default:            throw system::Exception ("Unknown mode in StorageFactory::createGroup");
    }
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline Partition<Type>* StorageFactory::createPartition (ICell* parent, const std::string& name, size_t nb)
{
    switch (_mode)
    {
        case STORAGE_FILE:  return StorageFileFactory::createPartition<Type> (parent, name, nb);
        case STORAGE_GZFILE:  return StorageGzFileFactory::createPartition<Type> (parent, name, nb);
        case STORAGE_COMPRESSED_FILE:  return StorageSortedFactory::createPartition<Type> (parent, name, nb);

        default:            throw system::Exception ("Unknown mode in StorageFactory::createPartition");
    }
}

/*********************************************************************
*********************************************************************/
template<typename Type>
inline CollectionNode<Type>* StorageFactory::createCollection (ICell* parent, const std::string& name, system::ISynchronizer* synchro)
{
    switch (_mode)
    {
        case STORAGE_FILE:  return StorageFileFactory::createCollection<Type> (parent, name, synchro);
        case STORAGE_GZFILE:  return StorageGzFileFactory::createCollection<Type> (parent, name, synchro);
        case STORAGE_COMPRESSED_FILE:  return StorageSortedFactory::createCollection<Type> (parent, name, synchro);

        default:            throw system::Exception ("Unknown mode in StorageFactory::createCollection");
    }
}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/
