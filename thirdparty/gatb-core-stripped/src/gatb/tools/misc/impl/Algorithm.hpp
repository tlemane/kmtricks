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

/** \file Algorithm.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Tool framework
 */

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_ALGORITHM_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_ALGORITHM_HPP_

/********************************************************************************/

#include <gatb/tools/designpattern/api/ICommand.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/api/IProperty.hpp>
#include <gatb/tools/misc/api/StringsRepository.hpp>
#include <gatb/tools/misc/impl/TimeInfo.hpp>
#include <gatb/tools/misc/impl/OptionsParser.hpp>
#include <gatb/tools/math/Integer.hpp>
#include <gatb/system/api/config.hpp>

#include <string>
#include <list>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/** \brief Framework class for implementing algorithm
 *
 * \see Tool
 */
class Algorithm : public system::SmartPointer
{
public:

    /** Constructor.
     * \param[in] name: name of the algorithm.
     * \param[in] nbCores : number of cores to be used for this algorithm.
     * \param[in] input : extra options for configuring the algorithm. */
    Algorithm (const std::string& name, int nbCores=-1, gatb::core::tools::misc::IProperties* input=0);

    /** Destructor. */
    virtual ~Algorithm ();

    /** Get tool name
     * \return the algorithm name. */
    std::string getName () const  { return _name; }

    /** Run the algorithm, ie. call 'execute'. */
    void run ();

    /** Execution of the algorithm. Abstract method, must be refined in subclasses. */
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
    virtual dp::IDispatcher*        getDispatcher ()  { return _dispatcher; }

    /** Get a TimeInfo instance for the tool. This object can be used for gathering
     * execution times of some parts of the \ref execute method.
     * \return the time info instance.
     */
    virtual TimeInfo&               getTimeInfo   ()  { return _timeInfo;   }

    /** Get information about operating system resources used during the execution.
     * \return operating system information.
     */
    virtual IProperties*            getSystemInfo ()  { return _systemInfo; }

    /** Create an iterator for the given iterator. If the verbosity is enough, progress bar information
     * can be displayed.
     * \param[in] iter : object to be encapsulated by a potential progress information
     * \param[in] nbIterations : number of iterations to be done.
     * \param[in] message : message used if progress information has to be displayed
     * \param[in] listener : listener to be used; if null, a new one is created
     * \return the created iterator.
     */
    template<typename Item> 
    dp::Iterator<Item>* createIterator (
        dp::Iterator<Item>* iter,
        size_t nbIterations=0,
        const char* message=0,
        dp::IteratorListener* listener = 0
    )
    {
        if (nbIterations > 0 && message != 0)
        {
            //  We create some listener to be notified every 1000 iterations and attach it to the iterator.
            if (listener == 0)  { listener = createIteratorListener (nbIterations, message); }

            dp::impl::SubjectIterator<Item>* iterSubject = new dp::impl::SubjectIterator<Item> (iter, nbIterations/100);
            iterSubject->addObserver (listener);

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


    /********************************************************************************/

    /** Utility function for easily running a kmers based algorithms.
     * It ensures that the correct instance of the provided functor is launched,
     * according to the kmer size (known at runtime). */
    template <template<size_t> class Functor>
    static int mainloop (tools::misc::IOptionsParser* parser, int argc, char* argv[])
    {
        // We get a handle on the provided parser.
        LOCAL (parser);

        try {
            // We parse the user options.
            tools::misc::IProperties* options = parser->parse (argc, argv);

            // We apply the functor that calls the correct implementation of the functor
            // according to the kmer size value.
            tools::math::Integer::apply<Functor> (options->getInt (STR_KMER_SIZE), options);
        }

        catch (tools::misc::impl::OptionFailure& e)  {  return e.displayErrors (std::cerr);                         }
        catch (system::Exception& e)                 {  std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;  }

        return EXIT_SUCCESS;
    }

protected:

    /** Computes the uri from an uri (ie add a prefix if any). */
    std::string getUriByKey (const std::string& key)  { return getUri (getInput()->getStr(key)); }

    /** Computes the uri from an uri (ie add a prefix if any). */
    std::string getUri (const std::string& str)  { return getInput()->getStr(STR_PREFIX) + str; }

    /** Setters. */
    void setInput      (IProperties*            input)       { SP_SETATTR (input);      }
    void setOutput     (IProperties*            output)      { SP_SETATTR (output);     }
    void setInfo       (IProperties*            info)        { SP_SETATTR (info);       }
    void setSystemInfo (IProperties*            systemInfo)  { SP_SETATTR (systemInfo); }
    void setDispatcher (dp::IDispatcher*        dispatcher)  { SP_SETATTR (dispatcher); }

private:

    /** Name of the tool (set at construction). */
    std::string _name;

    IProperties* _input;

    IProperties* _output;

    IProperties* _info;

    IProperties* _systemInfo;

    dp::IDispatcher* _dispatcher;

    /** */
    TimeInfo _timeInfo;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_ALGORITHM_HPP_ */
