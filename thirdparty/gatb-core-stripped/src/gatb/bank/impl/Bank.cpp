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

#include <gatb/bank/impl/Bank.hpp>

#include <gatb/bank/impl/BankFasta.hpp>

#include <gatb/bank/impl/BankBinary.hpp>
#include <gatb/bank/impl/BankAlbum.hpp>

#include <gatb/system/api/Exception.hpp>

using namespace std;

#define DEBUG(a)  //printf a

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
Bank::Bank ()
{
    /** We register most known factories. */
    _registerFactory_ ("album",  new BankAlbumFactory(),  false);
    _registerFactory_ ("fasta",  new BankFastaFactory(),  false);
    _registerFactory_ ("binary", new BankBinaryFactory(), false);

    DEBUG (("Bank::Bank,  found %ld factories\n", _factories.size()));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Bank::~Bank ()
{
    for (list<Entry>::iterator it = _factories.begin(); it != _factories.end(); it++)
    {
        (it->factory)->forget ();
    }
    _factories.clear();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Bank::_registerFactory_ (const std::string& name, IBankFactory* instance, bool beginning)
{
    /** We look whether the factory is already registered. */
    IBankFactory* factory = _getFactory_ (name);

    DEBUG (("Bank::registerFactory : name='%s'  instance=%p  => factory=%p \n", name.c_str(), instance, factory));

    if (factory == 0)
    {
        if (beginning)  { _factories.push_front (Entry (name, instance));  }
        else            { _factories.push_back  (Entry (name, instance));  }
        instance->use();
    }
    else
    {
        throw system::Exception ("Bank factory '%s already registered", name.c_str());
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
bool Bank::_unregisterFactory_ (const std::string& name)
{
    for (list<Entry>::iterator it = _factories.begin(); it != _factories.end(); it++)
    {
        if (it->name == name)  { if (it->factory)  { it->factory->forget(); }  _factories.erase(it);  return true; }
    }
    return false;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IBankFactory* Bank::_getFactory_ (const std::string& name)
{
    for (list<Entry>::iterator it = _factories.begin(); it != _factories.end(); it++)
    {
        if (it->name == name)  { return it->factory; }
    }
    return 0;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IBank* Bank::_open_ (const std::string& uri)
{
    DEBUG (("Bank::open : %s  nbFactories=%ld \n", uri.c_str(), _factories.size()));

    IBank* result = 0;
    for (list<Entry>::iterator it = _factories.begin(); result==0 && it != _factories.end(); it++)
    {
        result = it->factory->createBank(uri);
        DEBUG (("   factory '%s' => result=%p \n", it->name.c_str(), result ));

        if (result) // if one of the factories produce a result, we can just stop, no need to try other factories.
            break;
    }

    if (result == 0) { throw system::Exception ("Unable to open bank '%s' (if it is a list of files, perhaps some of the files inside don't exist)", uri.c_str()); }

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
std::string Bank::_getType_ (const std::string& uri)
{
    string result = "unknown";

    /** We try to create the bank; if a bank is valid, then we have the factory name. */
    for (list<Entry>::iterator it = _factories.begin(); it != _factories.end(); it++)
    {
        IBank* bank = it->factory->createBank(uri);
        if (bank != 0)
        {
            result = it->name;
			if(!result.compare("fasta"))
			{
				//distinguish fasta and fastq
				tools::dp::Iterator<Sequence>* its = bank->iterator(); LOCAL(its);
				its->first();
				if(!its->isDone())
				{
					std::string qual = its->item().getQuality();
					if(!qual.empty())
					{
						result= "fastq";
					}
				}
			}
			
            delete bank;
            break;
        }
    }

    return result;
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

