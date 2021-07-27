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

/** \file CollectionFile.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Collection interface
 *
 *  This file holds interfaces related to the Collection interface
 */

#ifndef _GATB_CORE_TOOLS_STORAGE_IMPL_COLLECTION_FILE_HPP_
#define _GATB_CORE_TOOLS_STORAGE_IMPL_COLLECTION_FILE_HPP_

/********************************************************************************/

#include <gatb/tools/collections/api/Collection.hpp>
#include <gatb/tools/collections/impl/BagFile.hpp>
#include <gatb/tools/collections/impl/IteratorFile.hpp>
#include <gatb/tools/collections/impl/CollectionAbstract.hpp>
#include <gatb/system/impl/System.hpp>
#include <json/json.hpp>

#include <string>
#include <vector>
#include <fstream>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace storage   {
namespace impl      {
/********************************************************************************/

/** \brief Implementation of the Collection interface with a file.
 *
 * This implementation reads/writes Item objects in a file.
 */
template <class Item> class CollectionFile : public collections::impl::CollectionAbstract<Item>, public system::SmartPointer
{
public:

    /** Constructor. */
    CollectionFile (const std::string& filename, size_t cacheItemsNb=10000)
        : collections::impl::CollectionAbstract<Item> (
             new collections::impl::BagFile<Item>(filename),
             new collections::impl::IterableFile<Item>(filename, cacheItemsNb)
             /* Note (Rayan): this isn't very clean. Two files objectss are opened, one by BagFile (in write mode) and one in IterableFile (in read mode).
              * With Clang/OSX, turns out the IterableFile was created before BagFile, causing some troubles.
              * Also this is opening the file twice, not nice. Anyway until I think of a better system, it's kept as it is, and IterableFile does a small hack*/
          ),  _name(filename), _propertiesName(filename+".props")
    {}

    /** Destructor. */
    virtual ~CollectionFile() {}

    /** \copydoc tools::collections::Collection::remove */
    void remove ()  {  
        gatb::core::system::impl::System::file().remove (_name);  
        gatb::core::system::impl::System::file().remove (_propertiesName);  
    }


    /* R: some code duplication with GroupFile, but it's the same in the h d f 5 case. not the best design.
     * so, collections can hold properties, so can groups.. */

    /** \copydoc tools::collections::Collection::addProperty */
    void addProperty (const std::string& key, const std::string value)
    {
        //std::cout << "CollectionFile addProperty called, key=" << key << " value=" << value<< std::endl;
        std::ifstream myfile (_propertiesName);
        std::string data, line;
        if (myfile.is_open())
        {
            while ( getline (myfile,line) )
                data += line;
            myfile.close();
        }
        json::JSON j;
        if (data.size() > 0)
            j = json::LoadJson(data);
        // otherwise json is empty and we create it
        
        j[key] = value;

        std::string s = j.dump();
        std::ofstream myfile2;
        myfile2.open (_propertiesName);
        myfile2 << s;
        myfile2.close();
    }

    /** \copydoc tools::collections::Collection::getProperty */
    std::string getProperty (const std::string& key)
    {
        //std::cout << "CollectionFile getProperty called, key=" << key << std::endl;
        std::string result;

        std::ifstream myfile (_propertiesName);
        std::string data, line;
        if (myfile.is_open())
        {
            while ( getline (myfile,line) )
                data += line;
            myfile.close();
        }
        json::JSON j;
        if (data.size() > 0)
        {
            j = json::LoadJson(data);
            result = j[key].ToString();
        }

        return result;
    }

private:

    std::string _name;
    std::string _propertiesName;
};

/********************************************************************************/
/* Experimental (not documented). */
template <class Item> class CollectionGzFile : public collections::impl::CollectionAbstract<Item>, public system::SmartPointer
{
public:
    
    /** Constructor. */
    CollectionGzFile (const std::string& filename, size_t cacheItemsNb=10000)
    : collections::impl::CollectionAbstract<Item> (
                                                   new collections::impl::BagGzFile<Item>(filename),
                                                   new collections::impl::IterableGzFile<Item>(filename, cacheItemsNb)
                                                   ),  _name(filename)
    {}
    
    /** Destructor. */
    virtual ~CollectionGzFile() {}
    
    /** \copydoc tools::collections::Collection::remove */
    void remove ()  {  gatb::core::system::impl::System::file().remove (_name);  }
    
private:
    
    std::string _name;
};
  
/********************************************************************************/
/* Experimental (not documented). */
template <class Item> class CollectionCountFile : public collections::impl::CollectionAbstract<Item>, public system::SmartPointer
{
public:
    
    /** Constructor. */
    CollectionCountFile (const std::string& filename, size_t cacheItemsNb=10000)
    : collections::impl::CollectionAbstract<Item> (
                                                   new collections::impl::BagCountCompressedFile<Item>(filename),
                                                   new collections::impl::IterableCountCompressedFile<Item>(filename, cacheItemsNb)                                                    ),  _name(filename)
    {}
    
    /** Destructor. */
    virtual ~CollectionCountFile() {}
    
    /** \copydoc tools::collections::Collection::remove */
    void remove ()  {  gatb::core::system::impl::System::file().remove (_name);  }
    
private:
    
    std::string _name;
};
    
/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_STORAGE_IMPL_COLLECTION_FILE_HPP_ */
