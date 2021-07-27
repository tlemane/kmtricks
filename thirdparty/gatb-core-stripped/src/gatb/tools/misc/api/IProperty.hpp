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

/** \file IProperty.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Interface for properties, ie. list of tag values with hierarchical feature
 */

#ifndef _GATB_CORE_TOOLS_MISC_IPROPERTY_HPP_
#define _GATB_CORE_TOOLS_MISC_IPROPERTY_HPP_

#include <gatb/system/api/ISmartPointer.hpp>
#include <stdlib.h>
#include <stdarg.h>
#include <list>
#include <set>
#include <iostream>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
/********************************************************************************/

/** \brief  Definition of a property as a [key,value] entry with a given depth.
 *
 * Define what a property can be. This is an extension of the [key,value] concept
 * since we add a notion a depth to this couple, which makes possible to have a tree
 * vision of a simple list of [depth,key,value] entries.
 *
 * Such instances are managed by the IProperties class, that acts as a container of IProperty
 * instances.
 *
 *  \see IProperties
 */
class IProperty : public system::SmartPointer
{
public:
    /** Constructor.
     * \param[in] aDepth : depth of the [key,value]
     * \param[in] aKey   : the key
     * \param[in] aValue : the value
     */
    IProperty (size_t aDepth, const std::string& aKey, const std::string& aValue)
        : depth(aDepth), key(aKey), value(aValue)  {}

    /** Constructor.
     * \param[in] aKey   : the key
     * \param[in] aValue : the value
     */
    IProperty (const std::string& aKey="", const std::string& aValue="")
        : depth(0), key(aKey), value(aValue)  {}

    /** Depth of the property. 0 should mean root property. */
    size_t      depth;

    /** Key of the property as a string. */
    std::string key;

    /** Value of the property as a string. */
    std::string value;

    /** Returns the value of the property as a C++ string.
     * \return the value
     */
    const std::string&  getValue  ()  { return value;                }

    /** Returns the value of the property as an integer (it supposes that the string represents an integer).
     * \return the value
     */
    long                getInt    ()  { return atol (value.c_str()); }

    /** Returns the value of the property as a float (it supposes that the string represents an float).
     * \return the value
     */
    double              getDouble    ()  { return atof (value.c_str()); }

    /** Returns the value of the property as a C string.
     * \return the value
     */
    const char*         getString ()  { return value.c_str();        }
};

/** An alias here... In the future, we should replace IProperty by Property. */
typedef IProperty  Property;

/********************************************************************************/

/** Visitor for a IProperty instance.
 *
 *  Define a Design Pattern Visitor for the IProperty instance.
 *
 *  In this case, we have only the IProperty class to visit (not a true classes hierarchy like one can
 *  find more classically for the Visitor DP), but we add two methods, one called before the IProperty
 *  visit, and one called after.
 *
 *  This can be seen as a improved way to iterate the IProperty items of a IProperties instance.
 *
 *  It is defined as a SmartPointer for easing instance life cycle management.
 */
class IPropertiesVisitor : public system::SmartPointer
{
public:

    /** Called before the true visit of the IProperty instance. */
    virtual void visitBegin    () = 0;

    /** Visit of the IProperty instance.
     * \param[in] prop : the instance to be visited.
     */
    virtual void visitProperty (IProperty* prop) = 0;

    /** Called after the true visit of the IProperty instance. */
    virtual void visitEnd      () = 0;
};

/********************************************************************************/

/** \brief Container of IProperty instances with DP Visitor capability.
 *
 *  This interface merely defines a container of IProperty instances; it contains
 *  several 'add' methods for adding IProperty instances into the container.
 *
 *  It is possible to retrieve a specific IProperty instance given a key.
 *
 *  The main method is 'accept'; its purpose is to visit each contained IProperty instance.
 *  Note that the only way to iterate the whole IProperty set is to define its own IPropertiesVisitor
 *  class and make it accepted by the IProperties instance; the 'visitProperty' method should be then
 *  called for each IProperty instance.
 *
 *  It is defined as a SmartPointer for easing instance life cycle management.
 *
 *  \see IProperty
 *  \see IPropertiesVisitor
 */
class IProperties : public system::SmartPointer
{
public:

    /** Accept a visitor (should loop over all IProperty instances).
     * \param[in] visitor : visitor to be accepted
     */
    virtual void accept (IPropertiesVisitor* visitor) = 0;

    /** Add a IProperty instance given a depth, a key and a value provided in a printf way.
     * \param[in] depth  : depth of the property to be added
     * \param[in] aKey   : key of the property to be added
     * \param[in] format : define the format of the value of the property, the actual value being defined by the ellipsis
     * \return a IProperty instance is created and returned as result of the method.
     *
     */
    virtual IProperty* add (size_t depth, const std::string& aKey, const char* format=0, ...) = 0;

    /** Add a IProperty instance given a depth, a key and a value.
     * \param[in] depth  : depth of the property to be added
     * \param[in] aKey   : key of the property to be added
     * \param[in] aValue : value (as a string) of the property to be added
     * \return a IProperty instance is created and returned as result of the method.
     */
    virtual IProperty* add (size_t depth, const std::string& aKey, const std::string& aValue) = 0;

    /** Add all the IProperty instances contained in the provided IProperties instance. Note that a depth is provided
     *  and is added to the depth of each added IProperty instance.
     * \param[in] depth : depth to be added to each depth of added instances.
     * \param[in] prop  : instance holding IProperty instances to be added
     */
    virtual void       add (size_t depth, IProperties* prop) = 0;

    /** Add all the IProperty instances contained in the provided IProperties instance. Note that a depth is provided
     *  and is added to the depth of each added IProperty instance.
     * \param[in] depth : depth to be added to each depth of added instances.
     * \param[in] prop  : instance holding IProperty instances to be added
     */
    virtual void       add (size_t depth, const IProperties& prop) = 0;

    /** */
    virtual void add (IProperty* prop, va_list args) = 0;

    /** Merge the IProperty instances contained in the provided IProperties instance.
     * \param[in] prop  : instance holding IProperty instances to be added
     */
    virtual void  merge (IProperties* prop) = 0;

    /** Returns the IProperty instance given a key.
     * \param[in] key : the key
     * \return the IProperty instance if found, 0 otherwise.
     */
    virtual IProperty* operator[] (const std::string& key) = 0;

    /** Returns the IProperty instance given a key.
     * \param[in] key : the key
     * \return the IProperty instance if found, 0 otherwise.
     */
    virtual IProperty* get (const std::string& key) const = 0;

    /** Get the value of a property given its key.
     * \param[in] key : the key of the property
     * \return the value of the key as a string.
     */
    virtual std::string getStr    (const std::string& key) const  = 0;

    /** Get the value of a property given its key.
     * \param[in] key : the key of the property
     * \return the value of the key as an integer
     */
    virtual int64_t     getInt    (const std::string& key) const = 0;

    /** Get the value of a property given its key.
     * \param[in] key : the key of the property
     * \return the value of the key as a double.
     */
    virtual double      getDouble (const std::string& key) const = 0;

    /** Set the value of a property given its key.
     * \param[in] key : the key of the property
     * \param[in] value : value to be set.
     */
    virtual void setStr    (const std::string& key, const std::string& value) = 0;

    /** Set the value of a property given its key.
     * \param[in] key : the key of the property
     * \param[in] value : value to be set.
     */
    virtual void setInt    (const std::string& key, const int64_t& value) = 0;

    /** Set the value of a property given its key.
     * \param[in] key : the key of the property
     * \param[in] value : value to be set.
     */
    virtual void setDouble (const std::string& key, const double& value) = 0;

    /** Clone the instance
     * \return the cloned instance.
     */
    virtual IProperties* clone () = 0;

    /** Distribute arguments that are comma separated list.
     * \return the list of distributed IProperties instances.
     */
    virtual std::list<IProperties*> map (const char* separator) = 0;

    /** Get the known keys.
     * \return the set of keys
     */
    virtual std::set<std::string> getKeys () = 0;

    /** Move the item (given its key) to the front of the container.
     * \param[in] key : the key of the item to be moved.
     */
    virtual void setToFront (const std::string& key) = 0;

    /** Output the properties object through an output stream
     * \param[in] s : the output stream
     * \param[in] p : the properties object to output
     * \return the modified output stream
     */
    friend std::ostream & operator<<(std::ostream & s, const IProperties& p)  {  p.dump(s);  return s;  }

    /** Get the properties as an XML string
     * \return the XML string.
     */
    virtual std::string getXML () = 0;

    /** Fill a Properties instance from an XML stream.
     * \param[in] stream: the stream to be read (file, string...) */
    virtual void readXML (std::istream& stream) = 0;

protected:

    /** */
    virtual void dump (std::ostream& s) const  = 0;
};

/********************************************************************************/

#define PROP_END  ((IProperty*)0)

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IPROPERTY_HPP_ */
