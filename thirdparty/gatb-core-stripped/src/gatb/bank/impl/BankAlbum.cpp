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

#include <gatb/bank/impl/BankAlbum.hpp>
#include <gatb/bank/impl/Bank.hpp>
#include <gatb/bank/impl/BankFasta.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/impl/Tokenizer.hpp>

#include <fstream>

using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;
using namespace gatb::core::tools::misc;

#define DEBUG(a)    //printf a
#define VERBOSE(a)

/********************************************************************************/
namespace gatb {  namespace core {  namespace bank {  namespace impl {
/********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankAlbum::BankAlbum (const std::string& name, bool deleteIfExists) : _name(name)
{
    /** We get a handle on the file given its name. */
    system::IFile* file = getFile (name, deleteIfExists ? "w+" : NULL);

    /** We check that the provided name exists in filesystem. */
    if (file != 0)
    {
        char buffer[256];

        while (file->isEOF() == false)
        {
            /** We init the buffer. */
            *buffer = 0;

            /** We get the current line. */
            int len = file->gets (buffer, sizeof(buffer));

            if (len > 1)
            {
                /** We remove the end of line. */
                if (buffer[len-1] == '\n')  {  buffer[len-1] = 0; }

                string bankUri = buffer;

                /** We check whether it is a mere file name or there is also a directory name. */
                if (isOnlyFilename(buffer) == true)
                {
                    /** We add the path name of the album file. */
                    bankUri = System::file().getDirectory (name) + "/" + bankUri;
                }

                /** We add a new bank. */
                BankComposite::addBank (Bank::open(bankUri));

                /** We memorize the uri of this bank. */
                _banksUri.push_back (bankUri);
            }
        }

        delete file;
    }
    else
    {
        throw Exception ("Unable to use file '%s'", name.c_str());
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankAlbum::BankAlbum (const std::vector<std::string>& filenames)
{
    for (vector<string>::const_iterator it = filenames.begin(); it != filenames.end(); ++it)
    {
        /** We add a new bank. */
        BankComposite::addBank (Bank::open(*it));

        /** We memorize the uri of this bank. */
        _banksUri.push_back (*it);
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool BankAlbum::isAlbumValid (const std::string& uri)
{
    bool result = true;

    FILE* file = fopen (uri.c_str(), "r");
    if (file != 0)
    {
        char buffer[256];

        while (fgets(buffer, sizeof(buffer), file) != 0)
        {
            VERBOSE (("BankAlbum::isAlbumValid  BEFORE buffer='%s'\n", buffer));

            int len = strlen(buffer);

            /** We skip last characters. */
            for (len-- ; len >= 0; len--)  { if (isspace (buffer[len]) != 0)  { buffer[len]=0; } }

            VERBOSE (("BankAlbum::isAlbumValid  AFTER  buffer='%s'\n", buffer));

            if (strlen(buffer) > 0)
            {
                string bankUri = buffer;

                /** We check whether it is a mere file name or there is also a directory name. */
                if (isOnlyFilename(bankUri) == true)
                {
                    /** We add the path name of the album file. */
                    bankUri = System::file().getDirectory (uri) + "/" + bankUri;
                }

                /** We check whether the file exists. */
                bool exists = System::file().doesExist(bankUri);

                DEBUG (("BankAlbum::isAlbumValid  bankUri='%s' bankUri.size()=%d  exists=%d  result=%d \n",
                        bankUri.c_str(), bankUri.size(), exists, result
                ));

                if (!exists) { result=false; break; }
            }
        }
        fclose (file);
    }
    else
    {
        result = false;
    }

    return result;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IBank* BankAlbum::addBank (const std::string& bankUri)
{
    IBank* result = 0;

    DEBUG (("BankAlbum::add '%s'\n", bankUri.c_str() ));

    /** We add the uri into the album file. */
    system::IFile* file = getFile (_name);

    if (file != 0)
    {
        /** We write the uri in the file. */
        file->print ("%s\n", bankUri.c_str());

        /** We create a new bank. */
        result = Bank::open(bankUri);

        /** We put it into the album. */
        BankComposite::addBank (result);

        /** We memorize the uri of this bank. */
        _banksUri.push_back (bankUri);

        /** Cleanup. */
        delete file;
    }

    return result;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IBank* BankAlbum::addBank (const std::string& directory, const std::string& bankName, bool output_fastq, bool output_gz)
{
    IBank* result = 0;

    DEBUG (("BankAlbum::add '%s'\n", bankName.c_str() ));

    /** We add the uri into the album file. */
    system::IFile* file = getFile (_name, "a+");

    if (file != 0)
    {
        /** We build the bank uri. */
        string bankUri = directory + "/" + bankName;

        /** We write the uri in the file. */
        file->print ("%s\n", bankName.c_str());

        /** We create a new FASTA bank. */
        result = new BankFasta (bankUri, output_fastq, output_gz);

        /** We put it into the album. */
        BankComposite::addBank (result);

        /** We memorize the uri of this bank. */
        _banksUri.push_back (bankUri);

        /** Cleanup. */
        delete file;
    }

    return result;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
system::IFile* BankAlbum::getFile (const std::string& name, const char* mode)
{
    /** We check whether the file already exists or not. */
    if (mode==NULL) { mode = System::file().doesExist(name) == true ? "r+" : "w+"; }

    /** We create a handle on the file. */
    return System::file().newFile (name, mode);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool BankAlbum::isOnlyFilename (const std::string& path)
{
    if (path.empty())  { throw Exception ("Bad '%s' path in isOnlyFilename", path.c_str());  }

    /** It may not be bullet proof... */
    return (path[0]!='/' && path[0]!='.');
}


/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankAlbum::remove ()
{
    BankComposite::remove();
    System::file().remove (_name);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IBank* BankAlbumFactory::createBank (const std::string& uri)
{
    /** We check whether the uri is a "multiple" bank, i.e. a list (comma separated) of FASTA banks URIs. */
    tools::misc::impl::TokenizerIterator it (uri.c_str(), ",");

    /** We count the item number. */
    vector<string> names;
    for (it.first(); !it.isDone(); it.next())  { names.push_back(it.item()); }

    /** FIRST CASE : a comma separated list of banks uri. */
    if (names.size() > 1)
    {
        DEBUG (("BankAlbumFactory::createBank : count>1 (%d)\n", names.size()));

        vector<IBank*> banks;
        for (size_t i=0; i<names.size(); i++)
        {
            DEBUG (("   %s\n", names[i].c_str()));

            /** We create a vector with the 'unitary' banks. */
            banks.push_back (Bank::open (names[i]));
        }
        /** We return a composite bank. */
        return new BankComposite (banks);
    }

    /** SECOND CASE : an album file. */
    else if (names.size() == 1)
    {
        bool isAlbumValid = BankAlbum::isAlbumValid(uri);

        DEBUG (("BankAlbumFactory::createBank : count==1  isAlbumValid=%d\n", isAlbumValid));

        if (isAlbumValid == true)  {  return new BankAlbum (uri);  }
    }

    /** Nothing worked => return 0. */
    return 0;
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
