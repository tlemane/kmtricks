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

#include <gatb/tools/misc/impl/Tool.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/Progress.hpp>
#include <gatb/tools/misc/impl/LibraryInfo.hpp>
#include <gatb/tools/designpattern/impl/Command.hpp>

#define DEBUG(a)  //printf a

using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;

using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;

/********************************************************************************/
namespace gatb {  namespace core { namespace tools {  namespace misc {  namespace impl {
/********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Tool::Tool (const std::string& name) : userDisplayHelp(0), _helpTarget(0),userDisplayVersion(0), _versionTarget(0), _name(name), _input(0), _output(0), _info(0), _parser(0), _dispatcher(0)
{
    setOutput (new Properties());

    setInfo  (new Properties());
   // _info->add (0, _name);

    /** We create an options parser. */
    setParser (new OptionsParser(name));

    getParser()->push_back (new OptionOneParam (STR_NB_CORES,    "number of cores",      false, "0"  ));
    getParser()->push_back (new OptionOneParam (STR_VERBOSE,     "verbosity level",      false, "1"  ));
	getParser()->push_back (new OptionNoParam (STR_VERSION, "version", false));
	getParser()->push_back (new OptionNoParam (STR_HELP, "help", false));

	
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Tool::~Tool ()
{
    setInput      (0);
    setOutput     (0);
    setInfo       (0);
    setParser     (0);
    setDispatcher (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Tool::displayVersion(std::ostream& os){
	LibraryInfo::displayVersion (os);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IProperties* Tool::run (int argc, char* argv[])
{
    DEBUG (("Tool::run(argc,argv) => tool='%s'  \n", getName().c_str() ));
    try
    {
        /** We parse the user parameters. */
        IProperties* props = getParser()->parse (argc, argv);

        /** We run the tool. */
        return run (props);
    }
    catch (OptionFailure& e)
    {
        e.displayErrors (std::cout);
        return NULL;
    }
	catch (ExceptionHelp& h)
	{
		if(userDisplayHelp!=NULL)
		{
			this->userDisplayHelp(_helpTarget);
		}
		else
		{
			h.displayDefaultHelp (std::cout);
		}
		return NULL;
	}
	catch (ExceptionVersion& v)
	{
		if(userDisplayVersion!=NULL)
		{
			this->userDisplayVersion(_versionTarget);
		}
		return NULL;
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
IProperties* Tool::run (IProperties* input)
{
    /** We keep the input parameters. */
    setInput (input);

    if (getInput()->get(STR_VERSION) != 0)
    {
    	displayVersion(cout);
        return _output;
    }

    /** We define one dispatcher. */
    if (_input->getInt(STR_NB_CORES) == 1)
    {
        setDispatcher (new SerialDispatcher ());
    }
    else
    {
        setDispatcher (new Dispatcher (_input->getInt(STR_NB_CORES)) );
    }

    /** We may have some pre processing. */
    preExecute ();

    /** We execute the actual job. */
    {
        //TIME_INFO (_timeInfo, _name);
        execute ();
    }

    /** We may have some post processing. */
    postExecute ();

    /** We return the output properties. */
    return _output;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Tool::preExecute ()
{
    /** We add a potential config file to the input properties. */
    _input->add (1, new Properties (/*System::info().getHomeDirectory() + "/." + getName() */));

    /** We may have to add a default prefix for temporary files. */
//    if (_input->get(STR_PREFIX)==0)  { _input->add (1, STR_PREFIX, "tmp.");  }

    /** set nb cores to be actual number of free cores, if was 0. */
    if (_input->getInt(STR_NB_CORES)<=0)  { _input->setInt (STR_NB_CORES, System::info().getNbCores());  }

//    /** We add the input properties to the statistics result. */
//    _info->add (1, _input);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Tool::postExecute ()
{
//    /** We add the time properties to the output result. */
//    _info->add (1, _timeInfo.getProperties ("time"));
//
//    /** We add the output properties to the output result. */
//    _info->add (1, "output");
//    _info->add (2, _output);

    /** We may have to dump execution information into a stats file. */
//    if (_input->get(STR_STATS_XML) != 0)
//    {
//        XmlDumpPropertiesVisitor visit (_info->getStr (STR_STATS_XML));
//        _info->accept (&visit);
//    }

    /** We may have to dump execution information to stdout. */
    if (_input->get(STR_VERBOSE) && _input->getInt(STR_VERBOSE) > 0)
    {
        RawDumpPropertiesVisitor visit;
        _info->accept (&visit);
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
dp::IteratorListener* Tool::createIteratorListener (size_t nbIterations, const char* message)
{
    switch (getInput()->getInt(STR_VERBOSE))
    {
        case 0: default:    return new IteratorListener ();
        case 1:             return new ProgressTimer (nbIterations, message);
        case 2:             return new Progress      (nbIterations, message);
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
ToolComposite::ToolComposite (const std::string& name) : Tool(name)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
ToolComposite::~ToolComposite ()
{
    for (list<Tool*>::iterator it = _tools.begin(); it != _tools.end(); it++)
    {
        (*it)->forget ();
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
IProperties* ToolComposite::run (int argc, char* argv[])
{
    vector<IProperties*> inputs;

    /** We first parse the options for all tools. */
    for (list<Tool*>::iterator it = _tools.begin(); it != _tools.end(); it++)
    {
#if 0
        /** We get the parameters from the current parser. */
        IProperties* input = (*it)->getParser()->parse (argc, argv);

        /** We add the input into the vector that gather the tools inputs. */
        inputs.push_back (input);
#else

        try
        {
            /** We parse the user parameters. */
            (*it)->getParser()->parse (argc, argv);


			IProperties* input =  (*it)->getParser()->getProperties() ;
            /** We add the input into the vector that gather the tools inputs. */
            inputs.push_back (input);
        }
        catch (OptionFailure& e)
        {
			IProperties* input =  (*it)->getParser()->getProperties() ;

			/** We add the input into the vector that gather the tools inputs. */
			inputs.push_back (input);
			
//            e.getParser().displayErrors (stdout);
//            e.getParser().displayHelp   (stdout);
//            return NULL;
        }
#endif
    }

    IProperties* output = 0;
    size_t idx = 0;
    for (list<Tool*>::iterator it = _tools.begin(); it != _tools.end(); it++, idx++)
    {
        /** We get the parameters from the current inputs entry. */
        IProperties* input = inputs[idx];

        /** We may have to add the output of the previous tool to the input of the current tool.
         *  WARNING! The output of the previous tool should have a bigger priority than the
         *  user parameters of the current tool.
         */
        IProperties* actualInput = 0;
        if (output != 0)
        {
            actualInput = new Properties();
            actualInput->add (1, output);   // output of the previous tool
            actualInput->add (1, input);    // input  of the previous tool
        }
        else
        {
            actualInput = input;
        }

        /** We run the tool and get a reference on its output. */
        output = (*it)->run (actualInput);

        /** We add the current tool info to the global properties. */
        _info->add (1, (*it)->getInfo());
    }

    /** We return the output properties. */
    return _output;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void ToolComposite::add (Tool* tool)
{
    if (tool)
    {
        tool->use ();
        _tools.push_back(tool);
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
void ToolComposite::execute ()
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void ToolComposite::preExecute ()
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void ToolComposite::postExecute ()
{
    /** We may have to dump execution information into a stats file. */
//    if (_input->get(Tool::STR_STATS_XML) != 0)
//    {
//        XmlDumpPropertiesVisitor visit (_info->getStr (Tool::STR_STATS_XML), false);
//        _info->accept (&visit);
//    }
}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/
