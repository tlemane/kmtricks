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

/** \file Tool.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Tool framework
 */

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_TOOL_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_TOOL_HPP_

/********************************************************************************/

#include <gatb/tools/designpattern/api/ICommand.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/tools/misc/api/IProperty.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/OptionsParser.hpp>
#include <gatb/tools/misc/api/StringsRepository.hpp>

#include <string>
#include <list>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/** \brief Framework abstract class for implementing tools (ie. binary tools).
 *
 * This class provides facilities for:
 *  - parsing command line
 *  - dispatching work on several threads
 *  - getting execution times
 *  - gathering statistics information
 *
 * The Tool sub classes must implement the execute method; this is the place where
 * the actual job of the tool has to be done.
 *
 * If the tool is launch with the run method and the famous [argc,argv] couple,
 * the command line [argc,argv] is parsed through an options parser and the recognized
 * options are available through the getInput method.
 *
 * The option parser should be configured (ie. adding specific options) during the constructor
 * of the Tool sub class. By default, some default options are attached to the Tool options
 * parser; for instance :
 *      - "-nb-cores" is available and enables to set the number of cores that can be used by the tool.
 *      - "-verbose" is available and can be used for having progress bar while iterating iterators.
 *
 * A dispatcher is available with the getDispatcher method. This dispatcher may have been
 * configured (ie. set the number of available cores) during the constructor.
 *
 * Example:
 * \snippet ToyTool.cpp  snippet1
 *
 * \see Algorithm
 */
class Tool : public system::SmartPointer
{
public:

    /** Constructor.
     * \param[in] name: name of the tool. */
    Tool (const std::string& name);

    /** Destructor. */
    virtual ~Tool ();

    /** Get tool name
     * \return the tool name. */
    std::string getName () const  { return _name; }

    /** Run the tool with input parameters provided as a IProperties instance
     * \param[in] input : input parameters
     * \return the parsed options as a IProperties instance
     */
    virtual IProperties* run (IProperties* input);

    /** Run the tool with input parameters provided as a couple [argc,argv]
     * \param[in] argc : number of arguments
     * \param[in] argv : array of arguments
     * \return the parsed options as a IProperties instance
     */
    virtual IProperties* run (int argc, char* argv[]);

    /** Subclasses must implement this method; this is where the actual job of
     * the tool has to be done.
     */
    virtual void execute () = 0;

    /** Get the parsed options as a properties instance
     * \return the parsed options.
     */
    virtual IProperties*            getInput      ()  { return _input;      }

    /** Get output results as a properties instance
     * \return the output results
     */
    virtual IProperties*            getOutput     ()  { return _output;     }

    /** Get statistics information about the execution of the tool
     * \return the statistics
     */
    virtual IProperties*            getInfo       ()  { return _info;       }

    /** Get an option parser configured with recognized options for the tool
     * \return the options parser instance
     */
    virtual IOptionsParser*         getParser     ()  { return _parser;     }

    /** Get a dispatched that can be used for parallelization. The option "-nb-cores" can
     * be used, and thus the provided number is used for configuring the dispatcher.
     * \return the dispatcher for the tool
     */
    virtual dp::IDispatcher*        getDispatcher ()  { return _dispatcher; }

    /** Get a TimeInfo instance for the tool. This object can be used for gathering
     * execution times of some parts of the \ref execute method.
     * \return the time info instance.
     */
    virtual TimeInfo&               getTimeInfo   ()  { return _timeInfo;   }

    /** Create an iterator for the given iterable. If the verbosity is enough, progress bar information
     * can be displayed.
     * \param[in] iterable : object that creates the iterator.
     * \param[in] message : message used if progress information has to be displayed
     * \return the created iterator.
     */
    template<typename Item> dp::Iterator<Item>* createIterator (collections::Iterable<Item>& iterable, const char* message=0)
    {
        int64_t nbItems = (iterable.getNbItems() >= 0 ? iterable.getNbItems() : iterable.estimateNbItems());
        return createIterator (iterable.iterator(), nbItems, message);
    }

    /** Create an iterator for the given iterator. If the verbosity is enough, progress bar information
     * can be displayed.
     * \param[in] iter : object to be encapsulated by a potential progress information
     * \param[in] nbIterations : number of iterations to be done.
     * \param[in] message : message used if progress information has to be displayed
     * \return the created iterator.
     */
    template<typename Item> dp::Iterator<Item>* createIterator (dp::Iterator<Item>* iter, size_t nbIterations=0, const char* message=0)
    {
        if (nbIterations > 0 && message != 0)
        {
            //  We create some listener to be notified every 1000 iterations and attach it to the iterator.
            dp::impl::SubjectIterator<Item>* iterSubject = new dp::impl::SubjectIterator<Item> (iter, nbIterations/100);
            iterSubject->addObserver (createIteratorListener (nbIterations, message));

            /** We assign the used iterator to be the subject iterator. */
            iter = iterSubject;
        }

        /** We return the result. */
        return iter;
    }

    /** Creates an iterator listener according to the verbosity level.
     * \param[in] nbIterations : number of iterations to be done
     * \param[in] message : progression message
     * \return an iterator listener.
     */
    virtual dp::IteratorListener* createIteratorListener (size_t nbIterations, const char* message);

    /** Displays information about the GATB library
     * \param[in] os : output stream used for dumping library information
     */
    virtual void displayVersion(std::ostream& os);

	
	
	/* let the user provide its own function to display Help
	 */
	void setHelp(void (*user_Help)(void * target)) {     if(user_Help != NULL) userDisplayHelp = user_Help; }
	void setHelpTarget(void * helpTarget) { _helpTarget = helpTarget;}

	/* let the user provide its own function to display Version
	 */
	void setVersion(void (*user_Version)(void * target)) {     if(user_Version != NULL) userDisplayVersion = user_Version; }
	void setVersionTarget(void * versionTarget) { _versionTarget = versionTarget;}

	
protected:

    /** */
    virtual void preExecute  ();
    virtual void postExecute ();

    /** Computes the uri from an uri (ie add a prefix if any). */
    std::string getUriByKey (const std::string& key)  { return getUri (getInput()->getStr(key)); }

    /** Computes the uri from an uri (ie add a prefix if any). */
    std::string getUri (const std::string& str)  { return getInput()->getStr(STR_PREFIX) + str; }

    /** Setters. */
    void setInput      (IProperties*            input)       { SP_SETATTR (input);      }
    void setOutput     (IProperties*            output)      { SP_SETATTR (output);     }
    void setInfo       (IProperties*            info)        { SP_SETATTR (info);       }
    void setParser     (IOptionsParser*         parser)      { SP_SETATTR (parser);     }
    void setDispatcher (dp::IDispatcher*        dispatcher)  { SP_SETATTR (dispatcher); }

protected:

	
	//pointer to function to display help
	void (*userDisplayHelp)(void * target);
	void * _helpTarget;
	
	//pointer to function to display help
	void (*userDisplayVersion)(void * target);
	void * _versionTarget;

	
    /** Name of the tool (set at construction). */
    std::string _name;

    IProperties* _input;

    IProperties* _output;

    IProperties* _info;

    IOptionsParser* _parser;

    dp::IDispatcher* _dispatcher;

    /** */
    TimeInfo _timeInfo;

    friend class ToolComposite;
};

/********************************************************************************/

/* DEPRECATED. */
class ToolComposite : public Tool
{
public:

    /** Constructor.
     * \param[in] name: name of the tool. */
    ToolComposite (const std::string& name = "tool");

    /** */
    ~ToolComposite ();

    /** */
    IProperties* run (int argc, char* argv[]);

    /** */
    void add (Tool* tool);

private:

    std::list<Tool*> _tools;

    /** */
    void execute ();
    void preExecute  ();
    void postExecute ();
};

/********************************************************************************/

/* DEPRECATED. */
class ToolProxy : public Tool
{
public:

    /** */
    ToolProxy (Tool* ref) : Tool("proxy"), _ref (ref)  {}

    /** */
    virtual IOptionsParser* getParser ()  {  return _ref->getParser();  }

    /** */
    virtual IProperties* getInput  ()  { return _ref->getInput();     }
    virtual IProperties* getOutput ()  { return _ref->getOutput();    }
    virtual IProperties* getInfo   ()  { return _ref->getInfo();      }

    /** */
    virtual dp::IDispatcher*    getDispatcher ()   { return _ref->getDispatcher(); }

    /** */
    virtual TimeInfo&    getTimeInfo ()   { return _ref->getTimeInfo (); }

    /** */
    Tool* getRef ()  { return _ref; }

private:
    Tool* _ref;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_TOOL_HPP_ */
