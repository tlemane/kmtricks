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

/** \file Bank.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief User front end for opening genomic banks in a generic way
 */

#ifndef _GATB_CORE_BANK_IMPL_BANK_FACTORY_HPP_
#define _GATB_CORE_BANK_IMPL_BANK_FACTORY_HPP_

/********************************************************************************/

#include <gatb/bank/api/IBank.hpp>

#include <string>
#include <list>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace bank      {
namespace impl      {
/********************************************************************************/

/** \brief Front end for managing IBank objects.
 *
 * The Bank class can be used as a front end for banks management in GATB.
 *
 * Actually, its main purpose is to provide IBank instances from a given URI (likely
 * to be a FASTA file for instance).
 *
 * By using this class, clients can open genomic banks without knowing their actual
 * type: they only rely on the IBank interface. This means that developing tools
 * this way make them independent of the format of the input bank (at least for the
 * supported input format by GATB).
 *
 * Today, the following factories are registered:
 *  1) BankAlbumFactory
 *  2) BankFastaFactory
 *  3) BankBinaryFactory
 *
 * During a call to 'open', each factory is tried (in the order of registration)
 * until a correct IBank object is returned; if no valid IBank is found, an exception
 * is thrown.
 *
 * Example of use:
 * \snippet bank16.cpp  snippet16_bank
 *
 * NOTE : In case of a brand new bank creation, there is no corresponding method to 'open'
 * because one has to know exactly the format of the bank (it is not possible in general
 * to deduce the format only by analyzing the URI string and not the content).
 *
 * This class is a Singleton (private constructor) and its static public members are
 * accessed to this singleton.
 *
 * For developers adding new implementations of the IBank interface, they should register a
 * IBankFactory instance for their new class. Doing so make their new bank format available
 * through the Bank::open method.
 */
class Bank
{
public:

    /** Open a bank and get a IBank instance. Since the instance is created by a new
     * statement, the client has to release the IBank when it is no longer used (otherwise
     * a memory leak will happen).
     * \param[in] uri : uri of the bank.
     * \return the IBank instance. */
    static IBank* open (const std::string& uri)  { return singleton()._open_ (uri); }

    /** In case of a composite bank, return the number of sub banks.
     * \return number of sub banks. */
    static size_t getCompositionNb (const std::string& uri)  {  IBank* bank = open (uri);  LOCAL (bank);  return bank->getCompositionNb();  }

    /** Get the type of the bank as a string
     * \param[in] uri : uri of the bank.
     * \return the bank type as a string. */
    static std::string getType (const std::string& uri)  { return singleton()._getType_(uri); }

    /** Register a new factory, associated with a name.
     * Note that the position of the registered factory in the list of factories may be important; for instance,
     * if we want to add a custom fasta file factory, we should add it at the beginning in order to be sure that
     * our custom data format will be used by default.
     * \param[in] name : name of the factory
     * \param[in] instance : IBank factory
     * \param[in] beginning : if true add the factory in the first position of the factories list, if false in the last position. */
    static void registerFactory (const std::string& name, IBankFactory* instance, bool beginning)  { singleton()._registerFactory_ (name, instance, beginning); }

    /** Unregister a factory given its name.
     * \param[in] name : name of the factory to be unregistered.
     * \return true if the factory has been unregistered, false otherwise. */
    static bool unregisterFactory (const std::string& name)  { return singleton()._unregisterFactory_ (name); }

    /** Get a factory for a given name. */
    static IBankFactory* getFactory (const std::string& name)  { return singleton()._getFactory_(name); }

private:

    /** Private due to singleton method. */
    Bank ();

    /** Destructor. */
    ~Bank ();

    /** The order of registration is important, so we can't rely on a map and have to use
     * a list of Entry (equivalent to <key,value> of a map). */
    struct Entry
    {
        Entry (const std::string& name, IBankFactory* factory) : name(name), factory(factory) {}
        std::string   name;
        IBankFactory* factory;
    };

    std::list<Entry> _factories;

    /** Singleton instance. */
    static Bank& singleton()  { static Bank instance; return instance; }

    /** Wrapper for 'open' method. */
    IBank* _open_ (const std::string& uri);

    /** Wrapper for 'getType' method. */
    std::string _getType_ (const std::string& uri);

    /** Wrapper for 'registerFactory' method. */
    void _registerFactory_ (const std::string& name, IBankFactory* instance, bool beginning);

    /** Wrapper for 'unregisterFactory' method. */
    bool _unregisterFactory_ (const std::string& name);

    /** Wrapper for 'getFactory' method. */
    IBankFactory* _getFactory_ (const std::string& name);
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IMPL_BANK_FACTORY_HPP_ */
