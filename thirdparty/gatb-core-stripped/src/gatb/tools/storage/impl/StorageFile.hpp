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

/** \file StorageFile.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Collection interface
 *
 *  This file holds interfaces related to the Collection interface
 */

#ifndef _GATB_CORE_TOOLS_STORAGE_IMPL_STORAGE_FILE_HPP_
#define _GATB_CORE_TOOLS_STORAGE_IMPL_STORAGE_FILE_HPP_

/********************************************************************************/

#include <cassert>
#include <gatb/tools/storage/impl/CollectionFile.hpp>
#include <json/json.hpp>
#include <iostream>
#include <fstream>

#define DEBUG_STORAGE(a)  //printf a

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace storage   {
namespace impl      {
/********************************************************************************/

    class GroupFile : public Group
    {
        public:
            GroupFile (Storage* storage, ICell* parent, const std::string& name)
                : Group(storage->getFactory(),parent,name), filename("")
            {

                std::string prefix = storage->getName();
                folder = prefix;
                if (!system::impl::System::file().isFolderEndingWith(prefix,"_gatb"))
				    folder += "_gatb/";

				// create folder if it doesn't exist
				if(!system::impl::System::file().doesExistDirectory(folder)){
					int ok = system::impl::System::file().mkdir(folder, 0755);
					if(ok != 0){  
                        std::cout << "Error: can't create output directory (" << folder<< ")\n" << " debug, doesexist:" << system::impl::System::file().doesExistDirectory(folder);
                    std::cout << "created directory " << folder << std::endl; // doesn't seem to be ever printed
					}   
				}
				/** We may need to create the h d f 5 group. Empty name means root group, which is constructed by default. */
                if (name.empty() == false)
                {
                    filename = folder + parent->getFullId('.') + std::string(".") + name;
                }
                else
                {
                    filename = folder + parent->getFullId('.') + std::string(".") + "json" /* i chose this name arbitrarily*/;
                }

                std::ifstream myfile (filename);
                std::string data, line;
                if (myfile.is_open())
                {
                    while ( getline (myfile,line) )
                    {
                        data += line;
                    }
                    myfile.close();
                }
                //std::cout << "data:" << data<< std::endl;
                if (data.size() > 0)
                    j = json::LoadJson(data); // populate json array;

                //std::cout << "GroupFile initialized" << std::endl;
            }

            /** */
            ~GroupFile()
            {
                //std::cout << "groupfile destructor called, removing folder " << folder << std::endl;
                system::impl::System::file().rmdir(folder); // hack to remove the trashme folers. I'd have liked to make that call in remove() but for some reason remove() isn't called
            }

            /** */
            void addProperty (const std::string& key, const std::string value)
            {
                //std::cout << "GRoupFile addProperty called with: " << key << " / " << value << std::endl;
                j[key] = value;
                flushJson();               
            }

            /** */
            std::string getProperty (const std::string& key)
            {
                //std::cout << "GroupFile getProperty called with: " << key << std::endl;
                std::string result;

                if (j.hasKey(key)) {
                    return j[key].ToString();
                }
                return "";
            }

            /* this is unused - also the json api we're using doesn't support it - h d f 5 needed it for stupid reasons */
            void delProperty (const std::string& key)
            {
                j[key]="";
                flushJson();
            }

            void setProperty (const std::string& key, const std::string value)
            {
                //std::cout << "GRoupFile setProperty called with: " << key << " / " << value << std::endl;
                j[key] = value;
                flushJson(); 
            }

            void flushJson()
            {
                std::string s = j.dump();
                std::ofstream myfile;
                myfile.open (filename);
                myfile << s;
                myfile.close();
            }

            // for some reason this isn't called. let me know if it is
            void remove ()
            {
                std::cout << "GroupFile remove called" << std::endl;
            }

        private:
            json::JSON j;
            std::string filename;
            std::string folder;
    };


/** \brief Factory used for storage of kind STORAGE_FILE
 */
class StorageFileFactory
{
public:

    /** Create a Storage instance.
     * \param[in] name : name of the instance to be created
     * \param[in] deleteIfExist : if the storage exits in file system, delete it if true.
     * \param[in] autoRemove : auto delete the storage from file system during Storage destructor.
     * \return the created Storage instance
     */
    static Storage* createStorage (const std::string& name, bool deleteIfExist, bool autoRemove)
    {
        DEBUG_STORAGE (("StorageFileFactory::createStorage  name='%s'\n", name.c_str()));
        return new Storage (STORAGE_FILE, name, autoRemove);
    }

    /** Tells whether or not a Storage exists in file system given a name
     * \param[in] name : name of the storage to be checked
     * \return true if the storage exists in file system, false otherwise.
     */
    static bool exists (const std::string& name)
    {
        return false;
    }

    /** Create a Group instance and attach it to a cell in a storage.
     * \param[in] parent : parent of the group to be created
     * \param[in] name : name of the group to be created
     * \return the created Group instance.
     */
    static Group* createGroup (ICell* parent, const std::string& name)
    {
        DEBUG_STORAGE (("StorageFileFactory::createGroup  name='%s'\n", name.c_str()));

        ICell* root = ICell::getRoot (parent);
        Storage* storage = dynamic_cast<Storage*> (root);
        assert (storage != 0);

        return new GroupFile (storage, parent, name);
    }

    /** Create a Partition instance and attach it to a cell in a storage.
     * \param[in] parent : parent of the partition to be created
     * \param[in] name : name of the partition to be created
     * \param[in] nb : number of collections of the partition
     * \return the created Partition instance.
     */
    template<typename Type>
    static Partition<Type>* createPartition (ICell* parent, const std::string& name, size_t nb)
    {
        DEBUG_STORAGE (("StorageFileFactory::createPartition  name='%s'\n", name.c_str()));

        ICell* root = ICell::getRoot (parent);
        Storage* storage = dynamic_cast<Storage*> (root);
        assert (storage != 0);

        // get the partitions folder and prefix
        std::string storage_prefix = storage->getName();
        std::string file_folder = storage_prefix;
        if (!system::impl::System::file().isFolderEndingWith(storage_prefix,"_gatb"))
            file_folder += "_gatb/";

        std::string full_path = file_folder;
        std::string parent_base = parent->getFullId('.');
        std::string base_name = parent_base;
        if (parent_base.size() > 0)
            base_name += std::string(".");  // because gatb's getBaseName is stupid and cuts after the last dot
        base_name += name; 

        full_path += base_name; // but then base_name might have a suffix like ".1" for partitions

        //std::cout <<"name: " << name << " filename " << full_path << " prefix " << base_name<< std::endl;

        if (nb == 0)
        {   // if nb is 0, it means we're opening partitions and not creating them, thus we need to get the number of partitions.

           int nb_partitions=0;
           for (auto filename : system::impl::System::file().listdir(file_folder))
            {
                if (!filename.compare(0, base_name.size(), base_name)) // startswith
                {
                    nb_partitions++;
                }
			}
            nb = nb_partitions;
            if (nb == 0)
            {
                std::cout << "error: could not get number of partition for " << name << " using StorageFile" << std::endl;
                exit(1);
            }
            //std::cout << "got " << nb << " partitions" << std::endl;
		}
        else
        {
            // else, if nb is set, means we're creating some partitions. let's delete all the previous ones to avoid wrongly counting 
            for (auto filename : system::impl::System::file().listdir(file_folder))
            {
                //std::cout <<"name: " << name << " comparing " << filename << " with prefix " << base_name << std::endl;
                if (!filename.compare(0, base_name.size(), base_name)) // startswith
                {
                    // some additional guard:
                    if (filename == "." ||filename == "..") continue;
                    system::impl::System::file().remove(file_folder + "/" + filename);
                    //std::cout << "deleting " << file_folder << "/" << filename << std::endl;
                }
            }
        }

        Partition<Type>* result = new Partition<Type> (storage->getFactory(), parent, name, nb);
        return result;
    }

    /** Create a Collection instance and attach it to a cell in a storage.
     * \param[in] parent : parent of the collection to be created
     * \param[in] name : name of the collection to be created
     * \param[in] synchro : synchronizer instance if needed
     * \return the created Collection instance.
     */
    template<typename Type>
    static CollectionNode<Type>* createCollection (ICell* parent, const std::string& name, system::ISynchronizer* synchro)
    {
        ICell* root = ICell::getRoot (parent);
        Storage* storage = dynamic_cast<Storage*> (root);
        assert (storage != 0);

        std::string storage_prefix = storage->getName();
        std::string folder = storage_prefix;
        if (!system::impl::System::file().isFolderEndingWith(storage_prefix,"_gatb"))
            folder += "_gatb/";

		// create folder if it doesn't exist
		if(!system::impl::System::file().doesExistDirectory(folder)){
			int ok = system::impl::System::file().mkdir(folder, 0755);
			if(ok != 0){  
                std::cout << "Error: can't create output directory (" << folder<< ")\n" << " debug, doesexist:" << system::impl::System::file().doesExistDirectory(folder);
			}   
                    std::cout << "created directory " << folder << std::endl;
		}

        /** We define the full qualified id of the current collection to be created. */
        std::string actualName = folder + parent->getFullId('.') + std::string(".") + name;

		DEBUG_STORAGE (("StorageFileFactory::createCollection  name='%s'  actualName='%s' \n", name.c_str(), actualName.c_str() ));

        return new CollectionNode<Type> (storage->getFactory(), parent, name, new CollectionFile<Type>(actualName));
    }
};

/********************************************************************************/
/* Experimental (not documented). */
class StorageGzFileFactory
{
public:
    
    /** */
    static Storage* createStorage (const std::string& name, bool deleteIfExist, bool autoRemove)
    {
        return new Storage (STORAGE_GZFILE, name, autoRemove);
    }
    
    /** */
    static bool exists (const std::string& name)
    {
        return false;
    }

    /** */
    static Group* createGroup (ICell* parent, const std::string& name)
    {
        ICell* root = ICell::getRoot (parent);
        Storage* storage = dynamic_cast<Storage*> (root);
        assert (storage != 0);
        
        return new Group (storage->getFactory(), parent, name);
    }
    
    /** */
    template<typename Type>
    static Partition<Type>* createPartition (ICell* parent, const std::string& name, size_t nb)
    {
        ICell* root = ICell::getRoot (parent);
        Storage* storage = dynamic_cast<Storage*> (root);
        assert (storage != 0);
        
        Partition<Type>* result = new Partition<Type> (storage->getFactory(), parent, name, nb);
        return result;
    }
    
    /** */
    template<typename Type>
    static CollectionNode<Type>* createCollection (ICell* parent, const std::string& name, system::ISynchronizer* synchro)
    {
        /** We define the full qualified id of the current collection to be created. */
        std::string actualName = std::string("tmp.") + name;
        
        ICell* root = ICell::getRoot (parent);
        Storage* storage = dynamic_cast<Storage*> (root);
        assert (storage != 0);
        
        return new CollectionNode<Type> (storage->getFactory(), parent, name, new CollectionGzFile<Type>(actualName));
    }
};

/********************************************************************************/
/* Experimental (not documented). */
class StorageSortedFactory
{
public:
    
    /** */
    static Storage* createStorage (const std::string& name, bool deleteIfExist, bool autoRemove)
    {
        return new Storage (STORAGE_COMPRESSED_FILE, name, autoRemove);
    }
    
    /** */
    static bool exists (const std::string& name)
    {
        return false;
    }

    /** */
    static Group* createGroup (ICell* parent, const std::string& name)
    {
        ICell* root = ICell::getRoot (parent);
        Storage* storage = dynamic_cast<Storage*> (root);
        assert (storage != 0);
        
        return new Group (storage->getFactory(), parent, name);
    }
    
    /** */
    template<typename Type>
    static Partition<Type>* createPartition (ICell* parent, const std::string& name, size_t nb)
    {
        ICell* root = ICell::getRoot (parent);
        Storage* storage = dynamic_cast<Storage*> (root);
        assert (storage != 0);
        
        Partition<Type>* result = new Partition<Type> (storage->getFactory(), parent, name, nb);
        return result;
    }
    
    /** */
    template<typename Type>
    static CollectionNode<Type>* createCollection (ICell* parent, const std::string& name, system::ISynchronizer* synchro)
    {
        /** We define the full qualified id of the current collection to be created. */
        std::string actualName = std::string("tmp.") + name;
        
        ICell* root = ICell::getRoot (parent);
        Storage* storage = dynamic_cast<Storage*> (root);
        assert (storage != 0);
        
        return new CollectionNode<Type> (storage->getFactory(), parent, name, new CollectionCountFile<Type>(actualName));
    }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_STORAGE_IMPL_STORAGE_FILE_HPP_ */
