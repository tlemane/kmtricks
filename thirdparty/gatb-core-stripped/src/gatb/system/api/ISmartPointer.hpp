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

/** \file ISmartPointer.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Smart Pointer Design Pattern interface
 *
 *  Define tools for easing life cycle of objects. Also known as smart pointer.
 */

#ifndef _GATB_CORE_SYSTEM_SMART_POINTER_HPP_
#define _GATB_CORE_SYSTEM_SMART_POINTER_HPP_

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace system    {
/********************************************************************************/

/** \brief Tool for managing instances life cycle
 *
 *  The goal of this class is to share easily objects between clients.
 *
 *  This class has an integer attribute that acts as a reference counter, i.e. a
 *  token counting how many clients are interested by the instance.
 *
 *  This is useful for sharing instances; if a client is interested by using an
 *  instance 'obj', she/he may call 'obj->use()' which will increase the internal
 *  token number. When the client is no more interested by using the instance 'obj',
 *  she/he may call 'obj->forget()', which will decrease the internal token.
 *
 *  When the token becomes 0, the instance is automatically destroyed.
 *
 *  Note that use() and forget() are virtual; it may happen (for singleton management
 *  for instance) that some subclass needs to refine them.
 *
 *  This pattern is often known as Smart Pointer.
 *
 *  Note that the STL provides its own smart pointer mechanism, known as auto_ptr.
 *  Here, our approach relies on subclassing instead of template use. The interest of our
 *  approach is to ease methods prototypes writing; with STL approach, one needs to uses
 *  every time auto_ptr<T> instead of only T, which can lower the readability.
 *
 *  On the other hand, our approach may be a little more dangerous than the STL approach
 *  since one has to be sure to forget an instance when needed. In our case, we use Smart
 *  Pointers mainly for attributes in class, so only have to be careful in constructors
 *  and destructors. Moreover, one can use the SP_SETATTR macro which eases this process.
 *  Note also the LOCAL macro that eases the local usage of an instance.
 *
 *  \see SP_SETATTR
 *  \see LOCAL
 */
class ISmartPointer
{
public:

    /** Destructor. */
    virtual ~ISmartPointer () {}

    /** Use an instance by taking a token on it */
    virtual void use () = 0;

    /** Forget an instance by releasing a token on it */
    virtual void forget () = 0;
};

/********************************************************************************/

/** \brief Implementation of the ISmartPointer interface
 *
 * This class implements also a security against concurrent access by several clients acting
 * from different threads. This is achieved by using intrinsics __sync_fetch_and_add and
 * __sync_fetch_and_sub in use and forget respectively.
 *
 * This class can't be instantiated since its default constructor is protected.
 *
 *  \see SP_SETATTR
 *  \see LOCAL
 */
class SmartPointer : public virtual ISmartPointer
{
public:
    /** \copydoc ISmartPointer::use */
    void use    ()
    {
        __sync_fetch_and_add (&_counterRef, 1 );
        // DON'T DO _counterRef++  because possible real time issues...
    }

    /** \copydoc ISmartPointer::forget */
    void forget ()
    {
        __sync_fetch_and_sub (&_counterRef, 1 );
        // DON'T DO _counterRef--  because possible real time issues...

        if (_counterRef<=0)  { delete this; }
    }

protected:

    /** Constructor. */
    SmartPointer () : _counterRef(0)  {}

    /** Destruction of the instance will be automatically called when the reference counter becomes null.
     *  The destructor is private in order to avoid to directly delete the instance without using the
     *  'use/forget' mechanism.
     */
    virtual ~SmartPointer ()  {}

private:
    /** Counts the number of references of the instance.*/
    int _counterRef;
};

/********************************************************************************/

/** \brief Local usage of SmartPointer instance.
 *
 * Small utility for locally getting a reference on a smart pointer. It creates
 *  a local (ie created in the execution stack) object that takes a reference on
 *  a smart pointer and get rid of it when the execution is out of the statements
 *  block holding the local object.
 *
 *  Note that the LocalObject class is very close to the std::auto_ptr. The first
 *  one uses inheritance, the second uses templates.
 *
 *  Sample:
 *  \code
 *  void foo ()
 *  {
 *     {
 *        // we create an instance of a class that inherits from SmartPointer and link it to a LocalObject
 *        LocalObject object (new MyClass ());
 *
 *        // we can access the referenced instance
 *        object.getPtr ();
 *     }
 *
 *     // Here, the MyClass instance should have been automatically deleted.
 *  }
 *  \endcode
 *
 *  \see LOCAL
 */
class LocalObject
{
public:
    /** Constructor.
     * \param[in] ptr : the instance we want locally manage.
     */
    LocalObject (ISmartPointer* ptr) : _ptr(ptr)  { if (_ptr)  {  _ptr->use();  } }

    /** Destructor. */
    ~LocalObject () { if (_ptr)  {  _ptr->forget ();  } }

    /** Getter on the referenced instance.
     * \return the referenced SmartPointer instance. */
    ISmartPointer* getPtr ()  { return _ptr; }

private:
    /** The SmartPointer instance we want local life cycle management. */
    ISmartPointer* _ptr;
};

/********************************************************************************/

/** Macro that creates an instance of type LocalObject whose name is the prefix '__' followed by the provided argument.
 *
 *  Sample:
 *  \code
 *  void foo ()
 *  {
 *     {
 *        // we create an instance of a class that inherits from SmartPointer
 *        MyClass* object = new MyClass ();
 *
 *        // we want that this object lives only inside the including statements block.
 *        // note that we use the LOCAL macro
 *        LOCAL (object);
 *     }
 *
 *     // Here, the object should have been automatically deleted.
 *  }
 *  \endcode
 *
 *  \see LocalObject
 */
#define LOCAL(object)  gatb::core::system::LocalObject __##object (object)

/** Macro that generates an instructions block that manages life cycle of a class attributes.
 *  It is dedicated to simply write clever setter methods.
 *
 * As a convention, the attribute name must begin by an underscore. Then the provided argument here
 * is the attribute name without this leading underscore.
 *
 * A typical use of this macro is the following:
 * \li in the constructor, use the default initialization of the attribute as 0
 * \li in the constructor body, use the smart setter with a provided argument
 * \li in the destructor body, use the smart setter with null argument.
 * \li in the private (or protected) part, defines a smart setter for the attribute by using the SP_SETATTR macro
 *
 * Sample of code:
 *  \code
 *  class MyClass
 *  {
 *  public:
 *       MyClass (SomeSmartPointerClass* ptr) : _ptr(0)  { setPtr (ptr); }
 *      ~MyClass ()  { setPtr (0); }
 *  private:
 *      SomeSmartPointerClass* _ptr;
 *      void setPtr (SomeSmartPointerClass* ptr)  { SP_SETATTR(ptr); }
 *  };
 *  \endcode
 *
 *  \see SmartPointer
 */
#define SP_SETATTR(a)  \
{  \
    if (_##a == a)  { return;  }   /* just to be sure that we don't reuse the same object */  \
    if (_##a != 0)  {  _##a->forget (); } \
    _##a = a; \
    if (_##a != 0)  {  _##a->use    (); } \
}

/********************************************************************************/

class SmartObject
{
public:

    SmartObject (ISmartPointer* ref=0) : _ref(0)  { setRef(ref); }

    virtual ~SmartObject ()  { setRef(0); }

    SmartObject (const SmartObject& o) : _ref(0)  { setRef(o._ref); }

    SmartObject& operator= (const SmartObject& o)  {  if (this != &o)  {  setRef (o._ref);  }  return *this; }

    void setRef (ISmartPointer* ref)  {  SP_SETATTR(ref); }
    ISmartPointer* getRef ()  {  return _ref; }

    bool hasRef () const { return _ref != 0; }

private:

    ISmartPointer* _ref;
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_SYSTEM_SMART_POINTER_HPP_ */
