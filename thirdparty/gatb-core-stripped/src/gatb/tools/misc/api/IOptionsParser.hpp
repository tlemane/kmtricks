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

/** \file IOptionsParser.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Tool for parsing command line arguments
 */

#ifndef _GATB_CORE_TOOLS_MISC_IOPTION_PARSER_HPP_
#define _GATB_CORE_TOOLS_MISC_IOPTION_PARSER_HPP_

/********************************************************************************/

#include <gatb/system/api/ISmartPointer.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/system/api/Exception.hpp>

#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
/********************************************************************************/

/** Forward declarations for the Visitor definition. */
namespace impl      {
    class OptionsParser;
    class Option;
}

/** Visitor design pattern for IOptionsParser. */
class IOptionsParserVisitor
{
public:
    /** Destructor. */
    virtual ~IOptionsParserVisitor() {}

    /** Visit a instance of OptionsParser
     * \param[in] object : the object to be visited
     * \param[in] depth : level of the visit */
    virtual void visitOptionsParser (impl::OptionsParser& object, size_t depth) = 0;

    /** Visit a instance of Option
     * \param[in] object : the object to be visited
     * \param[in] depth : level of the visit */
    virtual void visitOption (impl::Option& object, size_t depth) = 0;
};

/********************************************************************************/

/** \brief Parser interface that analyzes command line options.
 *
 * Client can use this class for registering command line options specifications
 * and then can use it for parsing some command line options, typically given
 * as arguments of the 'main' function.
 *
 * This interface is intended to be implemented as a Composite design pattern,
 * so we will have a 'leaf' implementation (see Option) and a 'composite'
 * implementation (see OptionsParser)
 */
class IOptionsParser : public system::SmartPointer
{
public:

    /*************************************************************/
    /*******************    General  methods   *******************/
    /*************************************************************/

    /** Get name. */
    virtual const std::string& getName () const = 0;

    /** Associate a name to the parser.
     * \param[in] name : the name of the parser. */
    virtual void setName (const std::string& name) = 0;

    /** Set visibility status. */
    virtual void setVisible (bool status) = 0;

    /** Get visibility status. */
    virtual bool isVisible() const = 0;

    /** Get help. */
    virtual const std::string& getHelp () const = 0;

    /** Set help
     * \param[in] help : the help string */
    virtual void setHelp (const std::string& help) = 0;

    /*************************************************************/
    /*******************    Parsing  methods   *******************/
    /*************************************************************/

    /** Perform the analyze of the arguments.
     * \param[in] argc : number of command line arguments.
     * \param[in] argv : table of arguments
     * \return object with information about the parsing.
     */
    virtual misc::IProperties* parse (int argc, char** argv) = 0;

    /** Perform the analyze of the arguments.
     * \param[in] s : string containing the options to be parsed
     * \return object with information about the parsing.
     */
    virtual misc::IProperties* parseString (const std::string& s) = 0;

    /** Return the properties found during parsing.
     * \return the parsed properties. */
    virtual misc::IProperties* getProperties ()  = 0;

    /** Tells whether an option has been seen during parsing.
     * \param[in] name : name of the option to be checked
     * \return true if option was found during parsing, false otherwise. */
    virtual bool saw (const std::string& name) const = 0;

    /*************************************************************/
    /*******************   Composite methods   *******************/
    /*************************************************************/

    /** Add a parser child at the back of known parsers.
     * \param[in] parser : the child parser
     * \param[in] expandDepth : while depth is less than expandDepth, put all the children and not the 'parser' itself.
     * \param[in] visibility : visibility status.
     */
    virtual void push_back (IOptionsParser* parser, size_t expandDepth=0, bool visibility=true) = 0;

    /** Add a parser child at the front of known parsers.
     * \param[in] parser : the child parser
     * \param[in] expandDepth : while depth is less than expandDepth, put all the children and not the 'parser' itself.
     * \param[in] visibility : visibility status.
     */
    virtual void push_front (IOptionsParser* parser, size_t expandDepth=0, bool visibility=true) = 0;

    /** Get a parser given its name.
     * \param[in] name : name of the parser to be retrieved
     * \return the parser instance if found, 0 otherwise. */
    virtual IOptionsParser* getParser (const std::string& name) = 0;

    /** Get the children parsers.
     * \return a list of parsers.*/
    virtual std::list<IOptionsParser*>& getParsers () = 0;

    /*************************************************************/
    /*********************   Miscellaneous   *********************/
    /*************************************************************/

    /** Return the default properties
     * \return the default properties. */
    virtual misc::IProperties* getDefaultProperties ()  = 0;

    /** Visitor design pattern. */
    virtual void accept (IOptionsParserVisitor& visitor, size_t depth=0) = 0;

    /** Structure that provides parsing results. */
    struct Result
    {
        /** Provides the properties found during parsing. */
        misc::impl::Properties properties;

        /** Provides the errors found during parsing. */
        std::list<std::string> errors;
    };
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IOPTION_PARSER_HPP_ */
