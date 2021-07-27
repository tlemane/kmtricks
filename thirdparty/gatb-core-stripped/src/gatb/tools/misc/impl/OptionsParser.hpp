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

/** \file OptionsParser.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Tool for parsing command line arguments
 */

#ifndef _GATB_CORE_TOOLS_MISC_OPTION_PARSER_HPP_
#define _GATB_CORE_TOOLS_MISC_OPTION_PARSER_HPP_

/********************************************************************************/

#include <gatb/tools/misc/api/IOptionsParser.hpp>
#include <gatb/tools/misc/api/StringsRepository.hpp>

#include <list>
#include <set>
#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/** \brief Implementation of the IOptionsParser interface.
 *
 * This implementation represents the 'composite' part of the Composite design pattern.
 *
 * Example:
 * \snippet optionsparser1.cpp  snippet1
 */
class OptionsParser : public IOptionsParser
{
public:

    /** Constructor. */
    OptionsParser (const std::string& name="", const std::string& help="");

    /** Destructor. */
    virtual ~OptionsParser ();

    /*************************************************************/
    /*******************    General  methods   *******************/
    /*************************************************************/

    /** \copydoc IOptionsParser::getName */
    const std::string& getName () const  { return _name; }

    /** \copydoc IOptionsParser::setName */
    void setName (const std::string& name)  { _name=name; }

    /** \copydoc IOptionsParser::setVisible */
    void setVisible (bool status)  { _visible=status; }

    /** \copydoc IOptionsParser::isVisible */
    bool isVisible() const  { return _visible; }

    /** \copydoc IOptionsParser::getHelp */
    const std::string& getHelp () const  { return _help; }

    /** \copydoc IOptionsParser::setHelp */
    void setHelp (const std::string& help)  { _help = help; }

    /*************************************************************/
    /*******************    Parsing  methods   *******************/
    /*************************************************************/

    /** \copydoc IOptionsParser::parse */
    misc::IProperties* parse (int argc, char** argv);

    /** \copydoc IOptionsParser::parseString */
    misc::IProperties* parseString (const std::string& s);

    /** \copydoc IOptionsParser::getProperties */
    misc::IProperties* getProperties ()  { return _properties; }

    /** \copydoc IOptionsParser::saw */
    bool saw (const std::string& name) const;

    /*************************************************************/
    /*******************   Composite methods   *******************/
    /*************************************************************/

    /** \copydoc IOptionsParser::push_back */
    void push_back (IOptionsParser* parser, size_t expandDepth=0, bool visibility=true);

    /** \copydoc IOptionsParser::push_front */
    void push_front (IOptionsParser* parser, size_t expandDepth=0, bool visibility=true);

    /** \copydoc IOptionsParser::getParser */
    IOptionsParser* getParser (const std::string& name);

    /** \copydoc IOptionsParser::getParsers */
    std::list<IOptionsParser*>& getParsers ()  { return _parsers; }

    /*************************************************************/
    /*********************   Miscellaneous   *********************/
    /*************************************************************/

    /** \copydoc IOptionsParser::getDefaultProperties */
    misc::IProperties* getDefaultProperties ();

    /** \copydoc IOptionsParser::accept */
    void accept (IOptionsParserVisitor& visitor, size_t depth=0)  { visitor.visitOptionsParser(*this, depth); }

protected:

    std::string _name;
    bool        _visible;
    std::string _help;

    std::list<IOptionsParser*> _parsers;

    misc::IProperties* _properties;
    void setProperties (misc::IProperties* properties)  { SP_SETATTR(properties); }

    std::ostream& indent (std::ostream& os, size_t level)  const { for (size_t i=0; i<level; i++)  { os << "   "; }  return os; }
};

/********************************************************************************/

/** \brief Implementation of the IOptionsParser interface.
 *
 * This implementation represents the 'leaf' part of the Composite design pattern.
 */
class Option : public OptionsParser
{
public:

    /** Constructor.
     * \param[in] name : name of the option
     * \param[in] nbArgs : number of arguments for this option
     * \param[in] mandatory : tells whether this option is mandatory or not
     * \param[in] defaultValue : default value for the option
     * \param[in] visible : tells whether this option may be shown in help
     * \param[in] help : textual help for this option
     */
    Option (
        const std::string&  name,
        int                 nbArgs,
        bool                mandatory,
        const std::string&  defaultValue,
        bool                visible,
        const std::string&  help
    )
        : OptionsParser(name, help), _nbArgs(nbArgs), _mandatory(mandatory), _defaultParam(defaultValue)
    {
        setVisible(visible);
    }

    /** Desctructor. */
    virtual ~Option() {}

    /** Get the default value for the option
     * \return the default value (may be empty) */
    std::string getDefaultValue () const { return _defaultParam; }

    /** Set the default value for the option
     * \param[in] value : the default value (may be empty) */
    void setDefaultValue (const std::string& value)  { _defaultParam = value; }

    /** Tells whether the option is mandatory or not.
     * \return the mandatory status.
     */
    bool isMandatory () const { return _mandatory; }

    /** \copydoc IOptionsParser::getParser */
    IOptionsParser* getParser (const std::string& name) { return name==getName() ? this : 0; }

    /** \copydoc IOptionsParser::accept */
    void accept (IOptionsParserVisitor& visitor, size_t depth=0)  { visitor.visitOption(*this, depth); }

protected:

    /** Gives the number of arguments that must follow the option.
     * \return the arguments number.
     */
    size_t getNbArgs () const   { return _nbArgs; }

    /** When an option is recognized in the arguments list, we look the number of waited args and put
     * them in a list of string objects. This is this list that is given as argument of the proceed() method
     * that mainly will affect the given args to the variable given to the instantiation of the
     * (derived class) Option.
     */
    virtual void proceed (const std::list<std::string>& args, IProperties& props) = 0;

    size_t          _nbArgs;
    bool            _mandatory;
    std::string     _defaultParam;

    friend struct ParserVisitor;
    friend class OptionsHelpVisitor;
};

/********************************************************************************/

/** \brief Option that has no argument.
 *
 * This is a special option (with no name) that memorize the arguments that are not
 * involved with a known option.
 */
class OptionNoParam : public Option
{
public:

    /** Constructor.
     * \param[in] name : name of the option
     * \param[in] help : textual help for this option
     * \param[in] mandatory : tells whether this option is mandatory or not
     * \param[in] visible : tells whether this option may be shown in help
     */
    OptionNoParam (
        const std::string&  name,
        const std::string&  help,
        bool                mandatory = false,
        bool                visible   = true
    )
        : Option (name, 0, mandatory, "", visible, help)
    {
    }

    /** \copydoc Option::proceed */
    void proceed (const std::list<std::string>& args, IProperties& props)
    {
        props.add (0, getName(), "");
    }
};

/********************************************************************************/

/** \brief Option that has one argument.
 *
 * This is a special option with only one argument.
 */
class OptionOneParam : public Option
{
public:

    /** Constructor.
     * \param[in] name : name of the option
     * \param[in] help : textual help for this option
     * \param[in] mandatory : tells whether this option is mandatory or not
     * \param[in] defaultValue : default value for the option
     * \param[in] visible : tells whether this option may be shown in help
     */
    OptionOneParam (
        const std::string&  name,
        const std::string&  help,
        bool                mandatory    = false,
        const std::string&  defaultValue = "",
        bool                visible      = true
    )
        : Option (name, 1, mandatory, defaultValue, visible, help)
    {
    }

    /** \copydoc Option::proceed */
    void proceed (const std::list<std::string>& args, IProperties& props)
    {
        props.add (0, getName(), args.front());
    }
};

/********************************************************************************/
/**************************           VISITORS         **************************/
/********************************************************************************/

/** \brief Visitor that iterates children of an OptionsParser instance. */
struct HierarchyParserVisitor : public IOptionsParserVisitor
{
    /** \copydoc IOptionsParserVisitor::visitOptionsParser */
    void visitOptionsParser (OptionsParser& object, size_t depth);
};

/** \brief Visitor that display help. */
class OptionsHelpVisitor : public IOptionsParserVisitor
{
public:

    /** Constructor.
     * \param[out] os : output stream where the output must be put. */
    OptionsHelpVisitor (std::ostream& os) : os(os),nameMaxLen(0) {}

    /** \copydoc IOptionsParserVisitor::visitOptionsParser */
    void visitOptionsParser (OptionsParser& object, size_t depth);

    /** \copydoc IOptionsParserVisitor::visitOption */
    void visitOption (Option& object, size_t depth);

private:

    std::ostream& os;
    size_t        nameMaxLen;
    std::ostream& indent (std::ostream& os, size_t level) const;
};

/** \brief Visitor that sets the visibility for a list of options. */
struct VisibilityOptionsVisitor : public IOptionsParserVisitor
{
    /** Constructor.
     * \param[in] visibility : status to be set.
     * \param[in] ... : list of labels of options to be modified; MUST BE terminated by a 0. */
    VisibilityOptionsVisitor (bool visibility, ...);

    /** \copydoc IOptionsParserVisitor::visitOptionsParser */
    void visitOptionsParser (OptionsParser& object, size_t depth);

    /** \copydoc IOptionsParserVisitor::visitOption */
    void visitOption (Option& object, size_t depth);

    bool                  _visibility;
    std::set<std::string> _names;
};

/********************************************************************************/

/** \brief Exception class to be used for option management error.
 *
 * This class should be thrown when something went wrong during options parsing.
 */
class OptionFailure
{
public:

    /** Constructor.
     * \param[in] parser : the parser that threw the exception.
     * \param[in] result : information gathered during the parsing. */
    OptionFailure (IOptionsParser* parser, IOptionsParser::Result result) :_parser(parser), _result(result)  {}

    /** Constructor.
     * \param[in] parser : the parser that threw the exception.
     * \param[in] msg : message to be displayed */
    OptionFailure (IOptionsParser* parser, const std::string& msg) :_parser(parser), _msg(msg)  {}

    /** Display information about the failure on the provided output stream
     * \param[out] os : output stream to be used
     * \return EXIT_FAILURE */
    int displayErrors (std::ostream& os) const;

protected:
    IOptionsParser*        _parser;
    IOptionsParser::Result _result;
    std::string            _msg;
};
	

	//specific exception to handle help and version
class ExceptionHelp
	{
	public:

		ExceptionHelp(IOptionsParser* parser) :_parser(parser) {}
		
		int displayDefaultHelp (std::ostream& os) const
		{
			OptionsHelpVisitor visitor (os);
			_parser->accept (visitor, 0);
			return EXIT_FAILURE;
		}
	protected:
		IOptionsParser*        _parser;
	};
	
class ExceptionVersion
	{
	};


/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_OPTION_PARSER_HPP_ */
