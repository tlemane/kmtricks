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

/** \file Property.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Progress feature
 */

#ifndef _GATB_CORE_TOOLS_MISC_IMPL_PROPERTY_HPP_
#define _GATB_CORE_TOOLS_MISC_IMPL_PROPERTY_HPP_

/********************************************************************************/

#include <gatb/tools/misc/api/IProperty.hpp>

#include <string>
#include <iostream>
#include <sstream>
#include <stack>
#include <cstdio>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/** \brief Implementation of IProperties interface.
 *
 *   This class provides a constructor that can read a file holding properties as [key,value] entries
 *   (one per line). This is useful for instance for managing a configuration file.
 *
 * Sample of use:
 * \code
 * void sample ()
 * {
 *      // we create a IProperties instance.
 *      IProperties* props = new Properties ();
 *
 *      // we add some entries. Note the different depth used: we have a root property having 3 children properties.
 *      props->add (0, "root", "");
 *      props->add (1, "loud",   "len=%d", 3);
 *      props->add (1, "louder", "great");
 *      props->add (1, "stop",   "[x,y]=[%f,%f]", 3.14, 2.71);
 *
 *      // we create some visitor for dumping the props into a XML file
 *      IPropertiesVisitor* v = new XmlDumpPropertiesVisitor ("/tmp/test.xml");
 *
 *      // we accept the visitor; after that, the output file should have been created.
 *      props->accept (v);
 * }
 * \endcode
 */
class Properties : public IProperties
{
public:

    /** Constructor. If a file path is provided, it tries to read [key,value] entries from this file.
     * \param[in] rootname : the file (if any) to be read
     */
    Properties (const std::string& rootname = "");

    /** Copy constructor. */
    Properties (const Properties& p);

    /** Destructor. */
    virtual ~Properties ();

    /** Affectation operator. */
    Properties& operator= (const Properties& p);

    /** \copydoc IProperties::accept */
    void accept (IPropertiesVisitor* visitor);

    /** \copydoc IProperties::add */
    IProperty* add (size_t depth, const std::string& aKey, const char* format=0, ...);

    /** \copydoc IProperties::add(size_t,const std::string&,const std::string&)  */
    IProperty* add (size_t depth, const std::string& aKey, const std::string& aValue);

    /** \copydoc IProperties::add(size_t,IProperties*)  */
    void add (size_t depth, IProperties* prop);

    /** \copydoc IProperties::add(size_t,IProperties*)  */
    void  add (size_t depth, const IProperties& prop);

    /**  */
    void add (IProperty* p, va_list args);

    /** \copydoc IProperties::merge  */
    void  merge (IProperties* prop);

    /** \copydoc IProperties::operator[]  */
    IProperty* operator[] (const std::string& key);

    /** \copydoc IProperties::get */
    IProperty* get (const std::string& key) const ;

    /** \copydoc IProperties::getStr */
    std::string getStr    (const std::string& key) const ;

    /** \copydoc IProperties::getInt */
    int64_t     getInt    (const std::string& key) const ;

    /** \copydoc IProperties::getDouble */
    double      getDouble (const std::string& key) const ;

    /** \copydoc IProperties::setStr */
    void setStr    (const std::string& key, const std::string& value);

    /** \copydoc IProperties::setInt */
    void setInt    (const std::string& key, const int64_t& value);

    /** \copydoc IProperties::setDouble */
    void setDouble (const std::string& key, const double& value);

    /** \copydoc IProperties::clone  */
    IProperties* clone ();

    /** \copydoc IProperties::map  */
    std::list<IProperties*> map (const char* separator);

    /** \copydoc IProperties::getKeys  */
    std::set<std::string> getKeys ();

    /** \copydoc IProperties::setToFront  */
    void setToFront (const std::string& key);

    /** Fill a Properties instance from an XML stream.
     * \param[in] stream: the stream to be read (file, string...) */
    void readXML (std::istream& stream);

    /** Fill a Properties instance from an XML string.
     * \param[in] xmlString: the XML string to be read (file, string...) */
    void readXML (const std::string& xmlString)  {  std::stringstream ss;  ss << xmlString; readXML (ss);  }

    /** \copydoc IProperties::getXML  */
    std::string getXML ();

    /* */
    void dump (std::ostream& s) const;

    /** Read a file holding [key,value] entries and add them through 'add' method.
     * \param[in] file : the file to be read
     */
    void readFile (const std::string& file);

private:

    /** List of IProperty instances. */
    std::list<IProperty*> _properties;

    /* */
    IProperty* getRecursive (const std::string& key, std::list<IProperty*>::const_iterator& it) const ;
};

/** Overload output stream operator. */
inline std::ostream& operator<< (std::ostream& os, const Properties& props)  { props.dump(os); return os; }

/********************************************************************************/

/* Factorization of common stuf for IPropertiesVisitor implementations. */
class AbstractOutputPropertiesVisitor : public IPropertiesVisitor
{
public:

     AbstractOutputPropertiesVisitor (std::ostream& aStream);
     AbstractOutputPropertiesVisitor (const std::string& filename);
    ~AbstractOutputPropertiesVisitor ();

protected:

    std::ostream* _stream;
    std::string   _filename;
};

/********************************************************************************/

/** \brief XML serialization of a IProperties instance.
 *
 *  This kind of visitor serializes into a file the content of a IProperties instance.
 *
 *  The output format is XML; the 'depth' attribute of each IProperty instance is used
 *  as a basis for building the XML tree.
 */
class XmlDumpPropertiesVisitor : public AbstractOutputPropertiesVisitor
{
public:

    /** Constructor.
     * \param[in] filename : uri of the file where to serialize the instance.
     * \param[in] propertiesAsRoot
     * \param[in] shouldIndent : tells whether we should use indentation
     */
    XmlDumpPropertiesVisitor (const std::string& filename, bool propertiesAsRoot=true, bool shouldIndent = true);

    /** Constructor.
     * \param[in] aStream : output stream
     * \param[in] propertiesAsRoot
     * \param[in] shouldIndent : tells whether we should use indentation
     */
    XmlDumpPropertiesVisitor (std::ostream& aStream, bool propertiesAsRoot=true, bool shouldIndent = true);

    /** Desctructor. */
    virtual ~XmlDumpPropertiesVisitor ();

    /** \copydoc IPropertiesVisitor::visitBegin */
    void visitBegin ();

    /** \copydoc IPropertiesVisitor::visitEnd */
    void visitEnd   ();

    /** \copydoc IPropertiesVisitor::visitProperty */
    void visitProperty (IProperty* prop);

private:

    /** The name of the serialization file. */
    std::string             _name;

    /** Some stack for XML production. */
    std::stack<std::string> _stack;

    /** */
    int _deltaDepth;

    /** Internals. */
    void pop    (size_t depth);

    /** Indentation of the XML file. */
    void indent (size_t n);

    bool _firstIndent;
    bool _shouldIndent;

    /** Method for writing into the file. */
    void safeprintf (const char* format, ...);
};

/********************************************************************************/

/** \brief Raw dump of a IProperties instance.
 *
 * This file format is simply a list of lines, with each line holding the key and the value
 * (separated by a space character). Note that the depth information is lost.
 *
 * This kind of output is perfect for keeping properties in a configuration file for instance.
 * This is used by a tool for its configuration file '.toolrc'
 */
class RawDumpPropertiesVisitor : public IPropertiesVisitor
{
public:

    /** Constructor.
     * \param os : output stream where the visitor can dump information. */
    RawDumpPropertiesVisitor (std::ostream& os = std::cout, int width=40, char sep=':');

    /** Desctructor. */
    virtual ~RawDumpPropertiesVisitor ();

    /** \copydoc IPropertiesVisitor::visitBegin */
    void visitBegin () {}

    /** \copydoc IPropertiesVisitor::visitEnd */
    void visitEnd   () {}

    /** \copydoc IPropertiesVisitor::visitProperty */
    void visitProperty (IProperty* prop);

private:
    std::ostream& _os;
    int _width;
    char _sep;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IMPL_PROPERTY_HPP_ */
