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

/** \file TimeInfo.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Tool for time measurement feature
 */

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_TIMETOOLS_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_TIMETOOLS_HPP_

/********************************************************************************/

#include <gatb/tools/misc/api/IProperty.hpp>
#include <gatb/system/api/ITime.hpp>
#include <gatb/system/impl/System.hpp>

#include <map>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/** \brief Tool for time statistics.
 *
 * This class provides methods for getting time information between two execution points.
 *
 * One can use a label for getting a specific time duration; it is possible later to get
 * this duration by giving the label.
 *
 * Example of use:
 * \code
 void foo ()
 {
     TimeInfo t;

     t.addEntry ("part1");
     // do something here
     t.stopEntry ("part1");

     t.addEntry ("part2");
     // do something here
     t.stopEntry ("part2");

     // now, we dump the duration of part1 and part2:
     cout << "part1: " << t.getEntryByKey("part1") << "  "
          << "part2: " << t.getEntryByKey("part2") << endl;
 }
 * \endcode
  */
class TimeInfo : public system::SmartPointer
{
public:

    /** Default constructor. */
    TimeInfo ();

	~TimeInfo();

    /** Constructor taking a time factory.
     * \param[in] aTime : the time factory to be used.
     */
    TimeInfo (system::ITime& aTime);

    /** Get the start time for a given label.
     * \param[in] name : the label
     */
    virtual void start (const char* name);

    /** Get the stop time for a given label.
     * \param[in] name : the label
     */
    virtual void stop (const char* name);

    /** Merge the content of the current time info with the provided one.
     * \param[in] ti : info to merged. */
    TimeInfo& operator+= (TimeInfo& ti);

    /** Didive the time (useful when gathered times from several threads).
     * \param[in] nb : divisor. */
    TimeInfo& operator/= (size_t nb );

    /** Provides (as a map) all got durations for each known label/
     * \return a map holding all retrieved timing information.
     */
    const std::map <std::string, u_int32_t>& getEntries ();

    /** Retrieve the duration for a given label.
     * \param[in] key : the label we want the duration for.
     * \return the duration.
     */
    u_int32_t getEntryByKey (const std::string& key);

    /** Retrieve the duration for a given label in seconds
     * \param[in] key : the label we want the duration for.
     * \return the duration.
     */
    double get (const std::string& key) { return (double)getEntryByKey(key) / 1000.0; }

    /** Creates and return as a IProperties instance the whole timing information.
     * \param[in] root : root name of the properties to be returned.
     * \return the created IProperties instance.
     */
    virtual tools::misc::IProperties* getProperties (const std::string& root);

private:

    system::ITime&  _time;
    std::map <std::string, u_int32_t>  _entriesT0;
    std::map <std::string, u_int32_t>  _entries;
	gatb::core::system::ISynchronizer* _synchro;
};

/********************************************************************************/

/** \brief Helper for time info statistics.
 *
 * This class allows to get the execution time within an instruction block.
 *
 * See also the TIME_INFO macro that eases its usage.
 *
 * Example:
 * \code
 void foo ()
 {
     TimeInfo t;

     {
         LocalTimeInfo local (t, "part1");

         // do something here
     }

     // now, we dump the exec time of the instruction block enclosing the LocalTimeInfo instance
     cout << "part1: " << t.getEntryByKey("part1") << "  " endl;
 }
 * \endcode
 *
 * */
class LocalTimeInfo
{
public:

    /** Constuctor
     * \param[in] ti : time info object to be used
     * \param[in] txt : key of the exec time to be got
     */
    LocalTimeInfo (TimeInfo& ti, const std::string& txt) : _ti(ti), _txt(txt)  {  _ti.start (_txt.c_str());  }

    /** Destructor. */
    ~LocalTimeInfo ()   {  _ti.stop (_txt.c_str());   }

private:
    TimeInfo&   _ti;
    std::string _txt;
};

#define  TIME_INFO(ti,txt)  gatb::core::tools::misc::impl::LocalTimeInfo TimeInfoTmp##__LINE__(ti,txt)

#define  TIME_START(ti,txt) gatb::core::tools::misc::impl::TimeInfo ti;  { ti.start (txt);
#define  TIME_STOP(ti,txt)   ti.stop(txt); }


/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_TIMETOOLS_HPP_ */
