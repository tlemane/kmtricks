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

#include <gatb/tools/misc/impl/Property.hpp>
#include <gatb/tools/misc/impl/XmlReader.hpp>
#include <gatb/tools/misc/impl/Tokenizer.hpp>
#include <gatb/system/impl/System.hpp>
#include <sstream>
#include <fstream>

#define DEBUG(a)  //printf a

using namespace std;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;
using namespace gatb::core::tools::dp;

#define PROP_BUFFER_SIZE (2*1024)

/********************************************************************************/
namespace gatb {  namespace core { namespace tools {  namespace misc {  namespace impl {
/********************************************************************************/

class InsertionVisitor : public IPropertiesVisitor
{
public:
    InsertionVisitor (size_t depth, IProperties* ref, set<string> keys) : _depth(depth), _ref(ref), _keys (keys) {}
    virtual ~InsertionVisitor() {}

    void visitBegin () {}
    void visitEnd   () {}

    void visitProperty (IProperty* prop)
    {
        if (_ref &&  prop)
        {
            /** We add the prop only if there is not already an existing identical prop. */
            if (_keys.find (prop->key) == _keys.end())
            {
                _ref->add ( prop->depth + _depth, prop->key, prop->value);
            }
        }
    }
private:
    size_t       _depth;
    IProperties* _ref;
    set<string>  _keys;
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Properties::Properties (const std::string& rootname)
{
    if (rootname.empty()==false)  { this->add (0, rootname); }
    // if (!initfile.empty())  {   readFile (initfile);  }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Properties::Properties (const Properties& p)
{
    for (std::list<IProperty*>::const_iterator it = p._properties.begin(); it != p._properties.end(); it++)
    {
        this->add ((*it)->depth, (*it)->key, (*it)->value);
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
Properties::~Properties ()
{
    for (std::list<IProperty*>::iterator it = _properties.begin(); it != _properties.end(); it++)
    {
        delete *it;
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
Properties& Properties::operator= (const Properties& p)
{
    if (this != &p)
    {
        for (std::list<IProperty*>::const_iterator it = p._properties.begin(); it != p._properties.end(); it++)
        {
            this->add ((*it)->depth, (*it)->key, (*it)->value);
        }
    }
    return *this;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IProperties* Properties::clone ()
{
    IProperties* result = new Properties ();

    for (std::list<IProperty*>::iterator it = _properties.begin(); it != _properties.end(); it++)
    {
        result->add ((*it)->depth, (*it)->key, (*it)->value);
    }

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
void Properties::accept (IPropertiesVisitor* visitor)
{
    visitor->visitBegin ();

    for (std::list<IProperty*>::iterator it = _properties.begin(); it != _properties.end(); it++)
    {
        visitor->visitProperty (*it);
    }

    visitor->visitEnd ();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IProperty* Properties::add (size_t depth, const std::string& aKey, const char* format, ...)
{
    IProperty* result = 0;

    if (format != 0)
    {
        char buffer[PROP_BUFFER_SIZE];
        va_list ap;
        va_start (ap, format);
        vsnprintf (buffer, sizeof(buffer), format, ap);
        va_end (ap);

        result = new IProperty (depth, aKey, buffer);
        _properties.push_back (result);
    }
    else
    {
        result = new IProperty (depth, aKey, "");
        _properties.push_back (result);
    }

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
IProperty* Properties::add (size_t depth, const std::string& aKey, const std::string& aValue)
{
    IProperty* result = new IProperty (depth, aKey, aValue);
    _properties.push_back (result);
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
void Properties::add (size_t depth, IProperties* properties)
{
    if (properties)
    {
        LOCAL (properties);

        /** We accept a visitor. */
        set<string>  nokeys;
        InsertionVisitor v (depth, this, nokeys);
        properties->accept (&v);
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
void Properties::add (size_t depth, const IProperties& properties)
{
    /** We accept a visitor. */
    set<string>  nokeys;
    InsertionVisitor v (depth, this, nokeys);
    ((IProperties&)properties).accept (&v);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Properties::add (IProperty* prop, va_list args)
{
    if (prop != 0)
    {
        LOCAL (prop);
        this->add (0, prop->key, prop->value);

        for (IProperty* p = 0;  (p = va_arg(args, IProperty*)) != 0; )
        {
            LOCAL (p);
            this->add (0, p->key, p->value);
        }
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
void Properties::merge (IProperties* properties)
{
    if (properties)
    {
        LOCAL (properties);

        /** We accept a visitor. */
        set<string>  nokeys;
        InsertionVisitor v (0, this, this->getKeys());
        properties->accept (&v);
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
IProperty* Properties::operator[] (const std::string& key)
{
    return get (key);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IProperty* Properties::get (const std::string& key) const
{
    IProperty* result = 0;

    list<IProperty*>::const_iterator it = _properties.begin();

    TokenizerIterator token (key.c_str(), ".");
    for (token.first(); !token.isDone(); token.next())
    {
        result = getRecursive (token.item(), it);
        if (!result) { break; }
    }

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
IProperty* Properties::getRecursive (const std::string& key, std::list<IProperty*>::const_iterator& it) const
{
    IProperty* result = 0;
    for (; !result  &&  it != _properties.end(); it++)
    {
        if (key.compare ((*it)->key)==0)    { result = *it; }
    }
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
std::string Properties::getStr (const std::string& key) const
{
    IProperty* prop = get (key);

    if (prop == 0)  {  throw Exception ("Empty property '%s'", key.c_str());  }

    return prop->getValue();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
int64_t Properties::getInt (const std::string& key) const
{
    IProperty* prop = get (key);

    if (prop == 0)  {  throw Exception ("Empty property '%s'", key.c_str());  }

    return prop->getInt();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
double Properties::getDouble (const std::string& key) const
{
    IProperty* prop = get (key);

    if (prop == 0)  {  throw Exception ("Empty property '%s'", key.c_str());  }

    return prop->getDouble();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Properties::setStr (const std::string& key, const std::string& value)
{
    IProperty* prop = get (key);

    if (prop == 0)  {  throw Exception ("Empty property '%s'", key.c_str());  }

    prop->value = value;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Properties::setInt (const std::string& key, const int64_t& value)
{
    IProperty* prop = get (key);

    if (prop == 0)  {  throw Exception ("Empty property '%s'", key.c_str());  }

    char buffer[64];
    snprintf (buffer, sizeof(buffer), "%ld", value);

    prop->value = buffer;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Properties::setDouble (const std::string& key, const double& value)
{
    IProperty* prop = get (key);

    if (prop == 0)  {  throw Exception ("Empty property '%s'", key.c_str());  }

    char buffer[64];
    snprintf (buffer, sizeof(buffer), "%f", value);

    prop->value = buffer;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Properties::readFile (const string& filename)
{
    /** We first check that the file exists. */
    if (System::file().doesExist(filename) == true)
    {
        IFile* file = System::file().newFile (filename, "r");
        if (file != 0)
        {
            char buffer[256];

            while (file->gets (buffer, sizeof(buffer) ) != 0)
            {
                char* key = strtok (buffer, " \t\n");
                if (key != 0  &&  isgraph(key[0]))
                {
                    char* value = key + strlen (key) + 1;

                    for ( ;  value; ++value)  {  if (*value != ' '  &&  *value != '\t')  { break; }  }

                    /** We remove the end of line. */
                    value[strlen(value)-1] = 0;

                    add (0, key, (value ? value : ""));
                }
            }

            delete file;
        }
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
void Properties::dump (std::ostream& s) const
{
    /** We dump some execution information. */
    RawDumpPropertiesVisitor visit (s);
    ((Properties*)this)->accept (&visit);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Properties::readXML (std::istream& stream)
{
    /** We create an XML reader. */
    XmlReader reader (stream);

    /** We create some specific observer class. */
    class XmlObserver : public IObserver
    {
    public:
        XmlObserver (Properties* ref) : _ref(ref), _depth(-1), _currentProperty(0)  {}

        void update (EventInfo* evt, ISubject* subject)
        {
            XmlTagOpenEvent* e1 = dynamic_cast<XmlTagOpenEvent*> (evt);
            if (e1)
            {
                /** We open a tag => increase the depth. */
                _depth ++;

                /** We add a property and keep a reference on it. */
                _currentProperty = _ref->add (_depth, e1->_name.c_str(), _txt.c_str());

                DEBUG (("XmlTagOpenEvent (%d): name=%s\n", _depth, e1->_name.c_str()));
                return;
            }

            XmlTagCloseEvent* e2 = dynamic_cast<XmlTagCloseEvent*> (evt);
            if (e2)
            {
                /** We close a tag => decrease the depth. */
                _depth --;

                DEBUG (("XmlTagCloseEvent(%d) : name=%s\n", _depth, e2->_name.c_str()));
                return;
            }

            XmlTagTextEvent* e3 = dynamic_cast<XmlTagTextEvent*> (evt);
            if (e3)
            {
                /** We find a text => set it to the current property.  */
                if (_currentProperty)
                {
                    _currentProperty->value = e3->_txt;

                    /** We don't want to bother with dummy text tag. */
                    _currentProperty = NULL;
                }

                DEBUG (("XmlTagTextEvent (%d) : txt=%s\n", _depth, e3->_txt.c_str()));
                return;
            }
        }

    private:
        Properties* _ref;
        string      _name;
        string      _txt;
        int         _depth;
        IProperty*  _currentProperty;
    };

    /** We remove all existing properties. */
    for (std::list<IProperty*>::iterator it = _properties.begin(); it != _properties.end(); it++)  {  delete *it;  }
    _properties.clear ();

    /** We attach this kind of observer to the reader. */
    XmlObserver observer (this);
    reader.addObserver (&observer);

    /** We read the stream. */
    reader.read ();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
std::list<IProperties*> Properties::map (const char* separator)
{
    list<IProperties*> result;
#if 0
    list <Iterator<char*>* > itList;

    for (list<IProperty*>::iterator it = _properties.begin(); it != _properties.end(); it++)
    {
        itList.push_back (new TokenizerIterator ((*it)->getString(), separator));
    }
    DEBUG (("Properties::map : nbIt=%ld\n", itList.size() ));

    CartesianIterator <char*> p (itList);

    for (p.first(); ! p.isDone(); p.next())
    {
        list<char*>& current = p.currentItem();

        IProperties* dup = this->clone();
        result.push_back (dup);

        list<char*>::iterator      itStr;
        list<IProperty*>::iterator itProps;

        for (itStr=current.begin(), itProps=_properties.begin();
             itStr!=current.end() && itProps!=_properties.end();
             itStr++, itProps++
        )
        {
            if (*itStr != 0 && (*itProps)->value.compare(*itStr) != 0)
            {
                DEBUG (("key='%s'  current='%s'  new='%s'\n", (*itProps)->key.c_str(), (*itProps)->value.c_str(), (*itStr)));
                dup->getProperty ((*itProps)->key)->value = *itStr;
            }
        }
    }

    /** Some cleanup. */
    for (list<Iterator<char*>*>::iterator it = itList.begin(); it != itList.end(); it++)
    {
        delete *it;
    }
#endif
    /** We return the result. */
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
set<string> Properties::getKeys ()
{
    set<string> result;

    for (list<IProperty*>::iterator it = _properties.begin(); it != _properties.end(); it++)
    {
        result.insert (result.end(), (*it)->key);
    }

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
void Properties::setToFront (const std::string& key)
{
    for (list<IProperty*>::iterator it = _properties.begin(); it != _properties.end(); it++)
    {
        if (key.compare ((*it)->key)==0)
        {
            /** We move the found key to the beginning of the container. */
            _properties.splice (_properties.begin(), _properties, it);
            break;
        }
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
string Properties::getXML ()
{
    stringstream ss;
    XmlDumpPropertiesVisitor v(ss, false);
    this->accept (&v);
    return ss.str();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
AbstractOutputPropertiesVisitor::AbstractOutputPropertiesVisitor (std::ostream& aStream)
    : _stream(0)
{
    /** A stream is provided, we keep a reference on it. */
    _stream = &aStream;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
AbstractOutputPropertiesVisitor::AbstractOutputPropertiesVisitor (const std::string& filename)
    : _stream(0), _filename(filename)
{
    if (_filename.empty() == false)
    {
        /** We create a file. */
        _stream = new fstream (_filename.c_str(), ios::out);
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
AbstractOutputPropertiesVisitor::~AbstractOutputPropertiesVisitor ()
{
    if (_filename.empty() == false)
    {
        delete _stream;
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
XmlDumpPropertiesVisitor::XmlDumpPropertiesVisitor (
    const std::string& filename,
    bool propertiesAsRoot,
    bool shouldIndent
)
    : AbstractOutputPropertiesVisitor (filename),
      _name (propertiesAsRoot ? "properties" : ""), _deltaDepth(0), _firstIndent(true), _shouldIndent(shouldIndent)
{
    /** We add the initial tag. */
    if (_name.empty() == false)
    {
        indent (0);
        safeprintf ("<%s>", _name.c_str());
    }

    _deltaDepth = _name.empty() == false ? 0 : 1;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
XmlDumpPropertiesVisitor::XmlDumpPropertiesVisitor (
    std::ostream& aStream,
    bool propertiesAsRoot,
    bool shouldIndent
)
    : AbstractOutputPropertiesVisitor (aStream),
      _name (propertiesAsRoot ? "properties" : ""), _deltaDepth(0), _firstIndent(true), _shouldIndent(shouldIndent)
{
    /** We add the initial tag. */
    if (_name.empty() == false)
    {
        indent (0);
        safeprintf ("<%s>", _name.c_str());
    }

    _deltaDepth = _name.empty() == false ? 0 : 1;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
XmlDumpPropertiesVisitor::~XmlDumpPropertiesVisitor ()
{
    /** We add the final tag. */
    if (_name.empty() == false)
    {
        indent (0);
        safeprintf ("</%s>", _name.c_str());
    }

    safeprintf ("\n");
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void XmlDumpPropertiesVisitor::visitBegin ()
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
void XmlDumpPropertiesVisitor::visitEnd ()
{
    /** We dump the remaining tags. */
    pop (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void XmlDumpPropertiesVisitor::visitProperty (IProperty* prop)
{
    if (prop != 0)
    {
        size_t actualDepth = prop->depth + 1;

        DEBUG (("XmlDumpPropertiesVisitor::visitProperty:  actualDepth=%ld stack.size=%ld '%s' \n",
            actualDepth, _stack.size(), prop->key.c_str()
        ));

        if (actualDepth > _stack.size())
        {
            indent (actualDepth);
            safeprintf ("<%s>%s", prop->key.c_str(), prop->value.c_str());
            _stack.push (prop->key);
        }

        else if (actualDepth == _stack.size())
        {
            safeprintf ("</%s>", _stack.top().c_str());
            _stack.pop();

            indent (actualDepth);
            safeprintf ("<%s>%s",  prop->key.c_str(),  prop->value.c_str() );
            _stack.push (prop->key);
        }

        else
        {
            pop (actualDepth);

            indent (actualDepth);
            safeprintf ("<%s>%s",  prop->key.c_str(),  prop->value.c_str() );
            _stack.push (prop->key);
        }
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
void XmlDumpPropertiesVisitor::pop (size_t depth)
{
    safeprintf ("</%s>", _stack.top().c_str());
    _stack.pop();

    while (_stack.size() >= depth && !_stack.empty())
    {
        indent (_stack.size());
        safeprintf ("</%s>", _stack.top().c_str());
        _stack.pop();
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
void XmlDumpPropertiesVisitor::indent (size_t n)
{
    if (!_shouldIndent)  { return; }

    if (!_firstIndent)  {  safeprintf ("\n");  }

    for (size_t i=1; i<=(n-_deltaDepth); i++)  {  safeprintf ("   ");  }

    _firstIndent = false;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void XmlDumpPropertiesVisitor::safeprintf (const char* format, ...)
{
    /** A safe printf method that check that the output file is ok. */
    if (_stream != 0)
    {
        char buffer[PROP_BUFFER_SIZE];

        va_list ap;
        va_start (ap, format);
        vsnprintf (buffer, sizeof(buffer), format, ap);
        va_end (ap);

        (*_stream) << buffer;
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
RawDumpPropertiesVisitor::RawDumpPropertiesVisitor (std::ostream& os, int width, char sep)
    : _os(os), _width(width), _sep(sep)
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
RawDumpPropertiesVisitor::~RawDumpPropertiesVisitor ()
{
    _os.flush();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void RawDumpPropertiesVisitor::visitProperty (IProperty* prop)
{
    int width = _width;

    char buffer[1024];

    string indent;
    for (size_t i=0; i<prop->depth; i++)  { indent += "    "; }

    if (prop->getValue().empty() == false)
    {
        snprintf (buffer, sizeof(buffer), "%s%-*s %c %s\n", indent.c_str(), width, prop->key.c_str(), _sep, prop->value.c_str());
    }
    else
    {
        snprintf (buffer, sizeof(buffer), "%s%-*s\n", indent.c_str(), width, prop->key.c_str());
    }

    _os << buffer;
}

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/
