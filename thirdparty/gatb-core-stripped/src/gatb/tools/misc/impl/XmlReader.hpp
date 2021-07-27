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

/** \file XmlReader.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief XML parsing
 */

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_XML_READER_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_XML_READER_HPP_

/********************************************************************************/

#include <gatb/tools/designpattern/impl/Observer.hpp>
#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/* Root class for all XML events. Inherits from dp::EventInfo.
 */
class XmlEvent : public dp::EventInfo
{
public:
    /** Constructor. */
    XmlEvent () : dp::EventInfo(0) {}
};

/* XML event corresponding to a tag opening.
 */
class XmlTagOpenEvent : public XmlEvent
{
public:
    /** Constructor.
     * \param[in] name : name of the tag. */
    XmlTagOpenEvent (const std::string& name) : _name(name) {}

    /** Tag name*/
    std::string _name;
};

/* XML event corresponding to a tag closing.
 */
class XmlTagCloseEvent : public XmlEvent
{
public:
    /** Constructor.
     * \param[in] name : name of the tag. */
    XmlTagCloseEvent (const std::string& name) : _name(name) {}

    /** Tag name*/
    std::string _name;
};

/* XML event corresponding to a text
 */
class XmlTagTextEvent : public XmlEvent
{
public:
    /** Constructor.
     * \param[in] txt : found text. */
    XmlTagTextEvent (const std::string& txt) : _txt(txt) {}

    /** The found text. */
    std::string _txt;
};

/* XML event corresponding to an attribute
 */
class XmlTagAttributeEvent : public XmlEvent
{
public:
    /** Constructor.
     * \param[in] name : name of the attribute
     * \param[in] value : value of the attribute
     */
    XmlTagAttributeEvent (const std::string& name, const std::string& value)
        : _name(name), _value(value) {}

    /** Name of the attribute. */
    std::string _name;

    /** Value of the attribute. */
    std::string _value;
};

/********************************************************************************/

/** \brief Simple implementation of an XML (SAX) parser
 *
 * This implementation considers that the reader is mainly a Subject, and therefore
 * sends notification to potential listeners as its parsing goes on.
 *
 * Any attached observer is likely to receive instances of subclasses of XMLEvent, and
 * has to do something in reaction.
 *
 * \code
 * // We create some XML listener class
 * class XmlListener : public IObserver
 * {
 * public:
 *      void update (EventInfo* evt, ISubject* subject)
 *      {
 *          XmlTagOpenEvent* open = dynamic_cast<XmlTagOpenEvent*> (evt);
 *          if (open)  {  cout << "open '" << open->_name << "'" << endl;  return;  }
 *
 *          XmlTagCloseEvent* close = dynamic_cast<XmlTagCloseEvent*> (evt);
 *          if (close)  {  cout << "close '" << close->_name << "'" << endl;  return;  }
 *
 *          XmlTagTextEvent* text = dynamic_cast<XmlTagTextEvent*> (evt);
 *          if (text)  {  cout << "text '" << text->_txt << "'" << endl;  return;  }
 *
 *          XmlTagAttributeEvent* attribute = dynamic_cast<XmlTagAttributeEvent*> (evt);
 *          if (attribute)  {  cout << "attribute: name='" << attribute->_name << "'  value='" << attribute->_value << "'" << endl;  return;  }
 *     }
 * };
 *
 * void foo ()
 * {
 *      // We define some input stream holding an XML string
 *      stringstream is (stringstream::in | stringstream::out);
 *      is << "<properties><progression><exec_percentage>100</exec_percentage><nb_alignments>4257</nb_alignments></progression>properties>";
 *
 *      // We create a reader with our input stream.
 *      XmlReader reader (is);
 *
 *      // We instantiate one observer
 *       IObserver* listener = new XmlListener();
 *
 *      // We attach the listener to the reader
 *      reader.addObserver (listener);
 *
 *      // We read the XML stream
 *      reader.read();
 *
 *      // We detach the listener from the reader.
 *      reader.removeObserver (listener);
 *  }
 * \endcode
 *
 * \see XmlTagOpenEvent
 * \see XmlTagCloseEvent
 * \see XmlTagTextEvent
 * \see XmlTagAttributeEvent
 */
class XmlReader : public dp::impl::Subject
{
public:

    /** Constructor.
     * \param[in] is : the input stream containing the XML stream. */
    XmlReader (std::istream& is);

    /** Destructor. */
    virtual ~XmlReader();

    /** Parse the input stream and possibly sends notifications to potential observers.
     */
    void read ();

private:

    /** */
    std::istream& _is;

    /** */
    void processMainState (std::istream& is);
    void processTagState  (std::istream& is);

    /** */
    void normalizeText (std::string& s);

    /** */
    void replace (std::string& str, const std::string& oldStr, const std::string& newStr);
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_XML_READER_HPP_ */

