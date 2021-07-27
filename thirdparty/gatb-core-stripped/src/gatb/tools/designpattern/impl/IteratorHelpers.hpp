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

/** \file IteratorHelpers.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Some helper classes for iteration
 */

#ifndef _GATB_CORE_DP_ITERATOR_IMPL_ITERATOR_HELPERS_HPP_
#define _GATB_CORE_DP_ITERATOR_IMPL_ITERATOR_HELPERS_HPP_

#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <set>
#include <list>
#include <boost/variant.hpp>

/********************************************************************************/
namespace gatb  {
namespace core  {
namespace tools {
namespace dp    {
namespace impl  {
/********************************************************************************/

/** \brief Null implementation of the Iterator interface.
 *
 * This implementation merely iterates over nothing. It may be useful to use such
 * an instance when we have to provide an Iterator instance but with nothing to iterate.
 */
template <class Item> class NullIterator : public Iterator<Item>
{
public:

    /** \copydoc  Iterator::first */
    void first() {}

    /** \copydoc  Iterator::next */
    void next()  {}

    /** \copydoc  Iterator::isDone */
    bool isDone() { return true; }

    /** \copydoc  Iterator::item */
    Item& item ()  { static Item item;  return item; }

    /** Destructor. */
    ~NullIterator() {}
};

/********************************************************************************/

/** \brief Iterator over two iterators.
 *
 * We define a "product" iterator for two iterators, i.e. it will loop each possible
 * couple of the two provided iterators.
 *
 *  It is useful for having only one loop instead of two loops. Note however that it
 *  may still be more efficient to have two loops. The CartesianIterator is just here
 *  for easing the product iteration on small sets.
 *
 * NOTE: most likely it doesn't work in combination with Dispatcher, see PairedIterator for how to fix
 *
 * Example:
 * \snippet iterators2.cpp  snippet1
 */
template <class T1, class T2> class ProductIterator : public Iterator < std::pair<T1,T2> >
{
public:

    /** Constructor.
     * \param[in] it1 : first iterator.
     * \param[in] it2 : second iterator.
     */
    ProductIterator (Iterator<T1>& it1, Iterator<T2>& it2)  : _it1(it1), _it2(it2), _isDone(false)  {  first(); }

    /** Destructor. */
    virtual ~ProductIterator ()  {}

    /** \copydoc Iterator::first */
    void first()
    {
        /** We go to the beginning of the two iterators. */
        _it1.first();
        _it2.first();

        /** We use a specific attribute to keep track of the finish status of the iterator.
         *  This is merely an optimization in order not to call too often the "isDone" method
         *  on the two iterators.  */
        _isDone = false;

        /* to make it work with dispatcher the fix should be something like:
             if (!isDone)
             {
             std::pair<T1,T2> t;
             t.first = _it1->item();
             t.second = _it2->item();
             *(this->_item) = t;
            }
            */
    }

    /** \copydoc Iterator::next */
    void next()
    {
        /** We go to the next item of the second iterator. */
        _it2.next ();

        /** We check that it is not done. */
        if (_it2.isDone())
        {
            /** We go to the next item of the first iterator. */
            _it1.next ();

            if (! _it1.isDone() )
            {
                /** The first iterator is not done, we can reset the second to its beginning. */
                _it2.first ();
            }
            else
            {
                /** The first iterator is also done, the product iterator is therefore done. */
                _isDone = true;
            }
        }
        /* to make it work with dispatcher the fix should e something like:
             if (!isDone)
             {
             std::pair<T1,T2> t;
             t.first = _it1->item();
             t.second = _it2->item();
             *(this->_item) = t;
            }
            */
    }

    /** \copydoc Iterator::isDone */
    bool isDone() { return _isDone; }

    /** \copydoc Iterator::item */
    std::pair<T1,T2>& item ()
    {
        _current.first  = _it1.item();
        _current.second = _it2.item();

        return _current;
    }

private:

    /** First iterator. */
    Iterator<T1>& _it1;

    /** Second iterator. */
    Iterator<T2>& _it2;

    /** Current item in the iteration. */
    std::pair<T1,T2> _current;

    /** Tells whether the iteration is finished or not. */
    bool _isDone;
};

/********************************************************************************/

/** \brief Iterator over two iterators.
 *
 * We define a an iterator for two iterators, by iterating each of the two iterators
 * and providing a pair of the two currently iterated items.
 *
 *  Example:
 * \snippet iterators7.cpp  snippet1
 *
 */
template <class T1, class T2=T1>
class PairedIterator : public Iterator < std::pair<T1,T2> >
{
public:

    /** Constructor.
     * \param[in] it1 : first iterator.
     * \param[in] it2 : second iterator.
     */
    PairedIterator (Iterator<T1>* it1, Iterator<T2>* it2)  : _it1(0), _it2(0), _isDone(true)
    {
        setIt1 (it1);   setIt2 (it2);
    }

    /** Destructor. */
    virtual ~PairedIterator ()  {  setIt1 (0);   setIt2 (0);  }

    /** \copydoc Iterator::first */
    void first()
    {
        _it1->first();
        _it2->first();

        _isDone = _it1->isDone() || _it2->isDone();
        
        // due to the way Dispatcher works, we need to update _item, not just _current
        if (_isDone == false)  { 
            std::pair<T1,T2> t;
            t.first = _it1->item();
            t.second = _it2->item();
            *(this->_item) = t;
         }

    }

    /** \copydoc Iterator::next */
    void next()
    {
        _it1->next ();  _it2->next ();

        _isDone = _it1->isDone() || _it2->isDone();
        if (_isDone == false)  {  
            std::pair<T1,T2> t;
            t.first = _it1->item();
            t.second = _it2->item();
            *(this->_item) = t;
         }

    }

    /** \copydoc Iterator::isDone */
    bool isDone() { return _isDone;  }

    /** \copydoc Iterator::item */
    std::pair<T1,T2>& item () {
        //std::pair<T1,T2> t; // can't do a temporary object because apparently it'll be destructed, eh
        _current.first  = _it1->item();
        _current.second = _it2->item();
        *(this->_item) = _current; // not sure if this is essential but i'm keeping it

        return _current;
 }

private:

    /** First iterator. */
    Iterator<T1>* _it1;
    void setIt1 (Iterator<T1>* it1)  { SP_SETATTR(it1); }

    /** Second iterator. */
    Iterator<T2>* _it2;
    void setIt2 (Iterator<T2>* it2)  { SP_SETATTR(it2); }

    /** Finish status. */
    bool _isDone;

    std::pair<T1,T2> _current;

};

/********************************************************************************/

/** \brief Factorization of code for the subject part of the Observer pattern.
 */
class AbstractSubjectIterator
{
public:

    /** Constructor. */
    AbstractSubjectIterator () : _hasListeners(false), _isStarted(false) {}

    /** Destructor. */
    ~AbstractSubjectIterator ()
    {
        /** We remove all observers. */
        for (std::set<IteratorListener*>::iterator it = _listeners.begin(); it != _listeners.end(); it++)
        {
            (*it)->forget();
        }
    }

    /** Add an observer to the iterator. Such an observer is provided as a functor.
     * \param[in] f : functor to be subscribed to the iterator notifications.
     */
    void addObserver    (IteratorListener* f)
    {
        if (f != 0)
        {
            f->use ();
            _listeners.insert (f);
            _hasListeners=true;
        }
    }

    /** Remove an observer from the iterator. Such an observer is provided as a functor.
     * \param[in] f : functor to be unsubscribed from the iterator notifications.
     */
    void removeObserver (IteratorListener* f)
    {
        /** We look whether the given functor is already known. */
        std::set<IteratorListener*>::iterator lookup = _listeners.find (f);
        if (lookup != _listeners.end())
        {
            (*lookup)->forget();
            _listeners.erase (lookup);
            _hasListeners = _listeners.empty() == false;
        }
    }

    /** Set a message to the observers. */
    void setMessage (const std::string& message)
    {
        /** We remove all observers. */
        for (std::set<IteratorListener*>::iterator it = _listeners.begin(); it != _listeners.end(); it++)
        {
            (*it)->setMessage (message);
        }
    }

protected:

    /** Notify all the subscribed functors.
     * \param[in] current : number of currently iterated items during the iteration. */
    void notifyInc (u_int64_t current)
    {
        if (_isStarted == true)
        {
            /** We call each subscribing functor. */
            for (std::set<IteratorListener*>::iterator it = _listeners.begin(); it != _listeners.end(); it++)
            {
                /** Not very pretty syntax, but that works. */
               (*it)->inc (current);
            }
        }
    }

    /** Notify all the subscribed functors about the start of the iteration. */
    void notifyInit ()
    {
        if (_isStarted == false)
        {
            _isStarted = true;

            /** We call each subscribing functor. */
            for (std::set<IteratorListener*>::iterator it = _listeners.begin(); it != _listeners.end(); it++)
            {
                (*it)->init ();
            }
        }
    }

    /** Notify all the subscribed functors about the start of the iteration. */
    void notifyFinish ()
    {
        if (_isStarted == true)
        {
            _isStarted = false;

            /** We call each subscribing functor. */
            for (std::set<IteratorListener*>::iterator it = _listeners.begin(); it != _listeners.end(); it++)
            {
                (*it)->finish ();
            }
        }
    }

private:

    std::set<IteratorListener*> _listeners;
    bool                        _hasListeners;
    bool                        _isStarted;
};

/********************************************************************************/

/** \brief Iterator that notifies some listener during iteration.
 *
 * Implementation note: we have to keep reference (through pointers) on functors
 * because we want them to be notified and not a copy of them.
 *
 * Note also that we don't allow to have twice the same observer (we use a set as
 * observers container).
 *
 * \snippet iterators4.cpp  snippet1_SubjectIterator
 */
template <class Item> class SubjectIterator : public Iterator<Item>, public AbstractSubjectIterator
{
public:

    /** Constructor
     * \param[in] ref : the referred iterator
     * \param[in] modulo : notifies every 'modulo' time
     * \param[in] listener : default listener attached to this subject (default value is 0)
     */
    SubjectIterator (Iterator<Item>* ref, u_int32_t modulo, IteratorListener* listener=0)
        : _ref(0), _current(0), _modulo(modulo==0 ? 1 : modulo)
    {
        /** We set the reference. */
        setRef (ref);

        /** We may have a listener. */
        if (listener != 0)  { this->addObserver (listener); }
    }

    /** Destructor. */
    ~SubjectIterator ()
    {
        /** We release the reference. */
        setRef (0);
    }

    /* \copydoc Iterator::first */
    void first ()
    {
        notifyInit ();
        _current = 0;
        _ref->first ();
    }

    /* \copydoc Iterator::isDone */
    bool isDone ()
    {
        bool res = _ref->isDone();  if (res)  { notifyFinish(); }  return res;
    }

    /* \copydoc Iterator::next */
    void next ()
    {
        _ref->next ();
        if ((_current % _modulo) == 0)  { notifyInc (_current);  _current=0; }
        _current++;
    }

    /* \copydoc Iterator::item */
    Item& item ()  { return _ref->item(); }

    /* */
    void setItem (Item& current)  { _ref->setItem(current); }

    /* */
    void reset ()  { _ref->reset(); }
	
	/* GR : this func was missing, previously caused subject iterator to return false composition*/
	std::vector<Iterator<Item>*> getComposition()  { return _ref->getComposition(); }
	


private:

    Iterator<Item>* _ref;
    void setRef (Iterator<Item>* ref)  { SP_SETATTR(ref); }

    u_int64_t _current;
    u_int64_t _modulo;
};

/********************************************************************************/
/** \brief Iterator that gather two iterators in a single loop.
 *
 * This iterator is equivalent to two iterators with one outer loop and one inner loop.
 */
template <typename T1, typename T2, typename Update> class CompoundIterator : public Iterator<T2>
{
public:
    /** Constructor.
     * \param[in] it1 : iterator for the outer loop.
     * \param[in] it2 : iterator for the inner loop.
     * \param[in] update : functor for updating inner loop with new item from outer loop
     */
    CompoundIterator (Iterator<T1>&  it1, Iterator<T2>& it2, const Update& update)  : _it1(it1), _it2(it2), _update(update) {}

    /** \copydoc Iterator::first */
    void first ()
    {
        _it1.first();
        if (!_it1.isDone())
        {
            _update (&_it2, &_it1.item());
            _it2.first ();
        }
    }

    /** \copydoc Iterator::next */
    void next ()
    {
        _it2.next ();
        if (_it2.isDone())
        {
            _it1.next();
            if (!_it1.isDone())
            {
                _update (&_it2, &_it1.item());
                _it2.first ();
            }
        }
    }

    /** \copydoc Iterator::isDone */
    bool isDone() { return _it1.isDone(); }

    /** \copydoc Iterator::item */
    T2& item()  { return _it2.item(); }

private:
    Iterator<T1>&  _it1;
    Iterator<T2>&  _it2;
    Update         _update;
};

/********************************************************************************/
/** \brief Iterator that truncate the number of iterations if needed.
 *
 * This iterator iterates a referred iterator and will finish:
 *      - when the referred iterator is over
 *   or - when a limit number of iterations is reached.
 *
 *  Example:
 * \snippet iterators3.cpp  snippet1
 */
template <class Item> class TruncateIterator : public Iterator<Item>
{
public:

    /** Constructor.
     * \param[in] ref : the referred iterator
     * \param[in] limit : the maximal number of iterations.
     * \param[in] initRef : will call 'first' on the reference if true
     */
    TruncateIterator (Iterator<Item>& ref, u_int64_t limit, bool initRef=true)
        : _ref(ref), _limit(limit), _currentIdx(0), _initRef(initRef), _isDone(true)  {}

    /** \copydoc  Iterator::first */
    void first()
    {
        _currentIdx=0;
        if (_initRef)  { _ref.first(); }

        /** We check whether the iteration is finished or not. */
        _isDone = _ref.isDone() || _currentIdx >= _limit;

        /** IMPORTANT : we need to copy the referred item => we just can't rely on a simple
         * pointer (in case of usage of Dispatcher for instance where a buffer of items is kept
         * and could be wrong if the referred items doesn't exist any more when accessing the
         * buffer). */
        if (!_isDone)  { *(this->_item) = _ref.item(); }
    }

    /** \copydoc  Iterator::next */
    void next()
    {
        _currentIdx++;
        _ref.next();

        /** We check whether the iteration is finished or not. */
        _isDone = _ref.isDone() || _currentIdx >= _limit;

        /** IMPORTANT : we need to copy the referred item => we just can't rely on a simple
         * pointer (in case of usage of Dispatcher for instance where a buffer of items is kept
         * and could be wrong if the referred items doesn't exist any more when accessing the
         * buffer). */
        if (!_isDone)  { *(this->_item) = _ref.item(); }
    }

    /** \copydoc  Iterator::isDone */
    bool isDone()  {  return _isDone;  }

    /** \copydoc  Iterator::item */
    Item& item ()  {  return *(this->_item);  }

private:

    Iterator<Item>& _ref;
    u_int64_t       _limit;
    u_int64_t       _currentIdx;
    bool            _initRef;
    bool            _isDone;
};


/********************************************************************************/
/** \brief Iterator that can be cancelled at some point during iteration
 *
 * This iterator iterates a referred iterator and will finish:
 *      - when the referred iterator is over
 *   or - when the cancel member variable is set to true
 *
 */
template <class Item> class CancellableIterator : public Iterator<Item>
{
public:

    /** Constructor.
     * \param[in] ref : the referred iterator
     * \param[in] initRef : will call 'first' on the reference if true
     */
    CancellableIterator (Iterator<Item>& ref, bool initRef=true)
        : _cancel(true),  _ref(ref),  _initRef(initRef), _isDone(true) {}

    /** \copydoc  Iterator::first */
    void first()
    {
        _cancel = false;

        if (_initRef)  { _ref.first(); }

        /** We check whether the iteration is finished or not. */
        _isDone = _ref.isDone();

        /** IMPORTANT : we need to copy the referred item => we just can't rely on a simple
         * pointer (in case of usage of Dispatcher for instance where a buffer of items is kept
         * and could be wrong if the referred items doesn't exist any more when accessing the
         * buffer).
         */
        if (!_isDone)  { *(this->_item) = _ref.item(); }
    }

    /** \copydoc  Iterator::next */
    void next()
    {
        _ref.next();

        /** We check whether the iteration is finished or not. */
        _isDone = _ref.isDone() || _cancel;

        /** IMPORTANT : we need to copy the referred item => we just can't rely on a simple
         * pointer (in case of usage of Dispatcher for instance where a buffer of items is kept
         * and could be wrong if the referred items doesn't exist any more when accessing the
         * buffer). */
        if (!_isDone)  { *(this->_item) = _ref.item(); }
    }

    /** \copydoc  Iterator::isDone */
    bool isDone()  {  return _isDone || _cancel;  }

    /** \copydoc  Iterator::item */
    Item& item ()  {  return *(this->_item);  }

    bool        _cancel;

private:

    Iterator<Item>& _ref;
    bool            _initRef;
    bool            _isDone;
};

/********************************************************************************/
/** \brief Iterator that filters out some iterated items
 *
 * This iterator iterates a referred iterator and will filter out some items according
 * to a functor provided at construction.
 *
 * Example:
 * \snippet iterators6.cpp  snippet1
 */
template <class Item, typename Filter> class FilterIterator : public ISmartIterator<Item>
{
public:

    /** Constructor.
     * \param[in] ref : the referred iterator
     * \param[in] filter : the filter on items. Returns true if item is kept, false otherwise. */
    FilterIterator (Iterator<Item>* ref, Filter filter) : _ref(0), _filter(filter), _rank(0)  { setRef(ref); }

    /** Destructor. */
    ~FilterIterator ()  { setRef(0); }

    /** \copydoc  Iterator::first */
    void first() { _rank=0; _ref->first(); while (!isDone() && _filter(item())==false) { _ref->next(); }  }

    /** \copydoc  Iterator::next */
    void next()  { _rank++; _ref->next();  while (!isDone() && _filter(item())==false) { _ref->next(); }  }

    /** \copydoc  Iterator::isDone */
    bool isDone() { return _ref->isDone();  }

    /** \copydoc  Iterator::item */
    Item& item ()  { return _ref->item(); }

    /** \copydoc  Iterator::setItem */
    void setItem (Item& i)  { _ref->setItem(i); }

    u_int64_t size () const  { return 0; }
    u_int64_t rank () const  { return _rank; }

private:

    Iterator<Item>* _ref;
    void setRef (Iterator<Item>* ref)  { SP_SETATTR(ref); }

    Filter      _filter;
    u_int64_t  _rank;
};

/********************************************************************************/

template <class Item> class VectorIterator : public Iterator<Item>
{
public:
    /** */
    VectorIterator (const std::vector<Item>& items) : _items(items), _idx(0), _nb (items.size())   {}

    /** */
    VectorIterator () : _idx(0), _nb (0)  {}

    /** */
    virtual ~VectorIterator () {}

    /** \copydoc  Iterator::first */
    void first()  {  _idx = -1;  next ();  }

    /** \copydoc  Iterator::next */
    void next()  { ++_idx;  if (_idx < _nb ) { this->_item = &(_items[_idx]); }  }

    /** \copydoc  Iterator::isDone */
    bool isDone() {   return _idx >= _nb;  }

    /** \copydoc  Iterator::item */
    Item& item ()  { return *(this->_item); }

protected:
    std::vector<Item> _items;
    int32_t           _idx;
    int32_t           _nb;
};

/********************************************************************************/

template <class Item> class VectorIterator2 : public Iterator<Item>
{
public:

    VectorIterator2 (std::vector<Item>& items) : _items(items), _idx(0), _nb (items.size())  {}

    /** */
    virtual ~VectorIterator2 () {}

    /** \copydoc  Iterator::first */
    void first()  {  _idx = -1;  next ();  }

    /** \copydoc  Iterator::next */
    void next()  { ++_idx;  if (_idx < _nb ) { *(this->_item) = (_items[_idx]); }  }

    /** \copydoc  Iterator::isDone */
    bool isDone() {   return _idx >= _nb;  }

    /** \copydoc  Iterator::item */
    Item& item ()  { return *(this->_item); }

protected:
    std::vector<Item>& _items;
    int32_t            _idx;
    int32_t            _nb;
};

/********************************************************************************/

template <template <class> class IteratorType , typename T1, typename T2=T1, typename T3=T2, typename T4=T3>
class IteratorVariant : public IteratorType <boost::variant<T1,T2,T3,T4> >
{
private:

    typedef boost::variant < IteratorType<T1>, IteratorType<T2>, IteratorType<T3>, IteratorType<T4> > Type;

    Type var;

    struct fct_first  : public boost::static_visitor<>     {  template<typename T>  void operator() (T& a) const { a.first();          }};
    struct fct_next   : public boost::static_visitor<>     {  template<typename T>  void operator() (T& a) const { a.next ();          }};
    struct fct_isDone : public boost::static_visitor<bool> {  template<typename T>  bool operator() (T& a) const { return a.isDone();  }};

    struct fct_set : public boost::static_visitor<>
    {
        boost::variant<T1,T2,T3,T4>& val;
        fct_set (boost::variant<T1,T2,T3,T4>& val) : val(val) {}
        template<typename T>  void operator() (T& a) const { val = a.item();  }
    };

    bool _isDone;

public:

    virtual ~IteratorVariant() {}

    template<typename T>
    IteratorVariant& operator=(const T& t)  {  var = t;  return *this;  }

    /**  */
    void first()
    {
        boost::apply_visitor (fct_first(), var);
        _isDone = boost::apply_visitor (fct_isDone(), var);
        if (!_isDone)  {  boost::apply_visitor (fct_set(*(this->_item)), var);  }
    }

    /** */
    void next()
    {
        boost::apply_visitor (fct_next(), var);
        _isDone = boost::apply_visitor (fct_isDone(), var);
        if (!_isDone)  {  boost::apply_visitor (fct_set(*(this->_item)), var);  }
    }

    /** */
    bool isDone() { return _isDone; }

    /** */
    boost::variant<T1,T2,T3,T4>& item ()  {  return *(this->_item);  }

    /** Get a reference on the object to be configured as the currently iterated item.
     * \param[in] i : object to be referred. */
    virtual void setItem (boost::variant<T1,T2,T3,T4>& i)  {  this->_item = &i;  }
};

/********************************************************************************/

/** \brief TO BE DONE...
 */
template <class Container, typename Type> class STLIterator : public Iterator<Type>
{
public:
    /** Constructor.
     * \param[in]  l : the list to be iterated
     */
    STLIterator (const Container& l)  : _l(l),_isDone(true) {}

    /** Destructor (here because of virtual methods). */
    virtual ~STLIterator ()  {}

    /** */
    void first()
    {
        _iter = _l.begin();
        _isDone = _iter == _l.end ();
        if (!_isDone)  { * this->_item = *_iter; }
    }

    /** */
    void next()
    {
        _iter++;
        _isDone = _iter == _l.end ();
        if (!_isDone)  { * this->_item = *_iter; }
    }

    /** */
    bool isDone() { return _isDone; }

    /** */
    Type& item()  { return * this->_item; }

private:

    /** List to be iterated. */
    Container _l;

    /** STL iterator we are going to wrap. */
    typename Container::iterator _iter;

    bool _isDone;
};

/********************************************************************************/

/** \brief Iterator that loops over std::list
 *
 *  This class is a wrapper between the STL list implementation and our own Iterator
 *  abstraction.
 *
 *  Note that the class is still a template one since we can iterate on list of anything.
 *
 *  \code
 *  void foo ()
 *  {
 *      list<int> l;
 *      l.push_back (1);
 *      l.push_back (2);
 *      l.push_back (4);
 *
 *      ListIterator<int> it (l);
 *      for (it.first(); !it.isDone(); it.next())       {   cout << it.currentItem() << endl;   }
 *  }
 *  \endcode
 */
template <class Type> class ListIterator : public STLIterator<std::list<Type>, Type>
{
public:
    ListIterator (const std::list<Type>& l)  :  STLIterator<std::list<Type>, Type> (l)  {}
};

/********************************************************************************/

template <class Type> class VecIterator : public STLIterator<std::vector<Type>, Type>
{
public:
    VecIterator (const std::vector<Type>& l)  :  STLIterator<std::vector<Type>, Type> (l)  {}
};

/********************************************************************************/

/** \brief Composite iterator
 *
 * This iterator takes a list of iterators as input and iterates each one of these
 * iterators.
 *
 *  Example:
 * \snippet iterators9.cpp  snippet1
 */
template <class Item>
class CompositeIterator : public Iterator <Item>
{
public:

    /** Constructor.
     * \param[in] iterators : the iterators vector
     */
    CompositeIterator (std::vector <Iterator<Item>*>&  iterators)
        : _iterators(iterators), _currentIdx(0), _currentIt(0), _isDone(true)
    {
        for (size_t i=0; i<_iterators.size(); i++)  { _iterators[i]->use(); }

        _currentIt = _iterators[_currentIdx];
    }

    /** Destructor. */
    virtual ~CompositeIterator ()
    {
        for (size_t i=0; i<_iterators.size(); i++)  { _iterators[i]->forget(); }
    }

    /** \copydoc Iterator::first */
    void first()
    {
        /** We initialize attributes. */
        _currentIdx = 0;
        _isDone     = true;

        /** We look for the first non finished iterator. */
        update (true);
    }

    /** \copydoc Iterator::next */
    void next()
    {
        _currentIt->next();
        _isDone = _currentIt->isDone();

        if (_isDone == true)  {  update (false);  }
    }

    /** \copydoc Iterator::isDone */
    bool isDone() { return _isDone; }

    /** \copydoc Iterator::item */
    Item& item ()  {  return _currentIt->item(); }

    /** IMPORTANT : the Item argument provided to 'setItem' must be the object to be modified by
     * one of the delegate iterator AND NOT the current item of CompositeIterator. Therefore,
     * we make point the delegate current item to this provided Item argument. */
    void setItem (Item& i)  {  _currentIt->setItem (i);  }

    /** Get a vector holding the composite structure of the iterator. */
    virtual std::vector<Iterator<Item>*> getComposition() { return _iterators; }

private:

    std::vector <Iterator<Item>*>  _iterators;

    size_t          _currentIdx;
    Iterator<Item>* _currentIt;

    bool _isDone;

    void update (bool isFirst)
    {
        if (_currentIdx >= _iterators.size()) { _isDone=true;  return; }

        if (!isFirst)  { _currentIdx++; }

        while ((int)_currentIdx<(int)_iterators.size() && _isDone == true)
        {
            Iterator<Item>* previous = _currentIt;

            /** We get the next iterator. */
            _currentIt = _iterators[_currentIdx];
            assert (_currentIt != 0);

            /** We have to take the reference of the previous iterator. */
            _currentIt->setItem (previous->item());

            /** We have to "first" this iterator. */
            _currentIt->first();

            /** We update the 'isDone' status. */
            _isDone = _currentIt->isDone();

            if (_isDone==true) { _currentIdx++; }

            /** We can finish the previous item (only if not first call). */
            if (!isFirst) { previous->finalize(); }
        }
    }
};

	
/********************************************************************************/
/** \brief Iterator adaptation from one type to another one
 */
template <class T1, class T2, class Adaptor>
class IteratorAdaptor : public Iterator <T2>
{
public:

    /** Constructor. */
    IteratorAdaptor (Iterator<T1>* ref) : _ref(0) { setRef(ref); }

    /** Destructor. */
    ~IteratorAdaptor ()  { setRef(0); }

    /** Method that initializes the iteration. */
    void first() {  _ref->first(); }

    /** Method that goes to the next item in the iteration.
     * \return status of the iteration
     */
    void next()  { _ref->next(); }

    /** Method telling whether the iteration is finished or not.
     * \return true if iteration is finished, false otherwise.
     */
    bool isDone()  { return _ref->isDone(); }

    /** */
    T2& item ()  { return Adaptor() (_ref->item()); }

private:

    Iterator<T1>* _ref;
    void setRef (Iterator<T1>* ref)  { SP_SETATTR(ref); }
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_DP_ITERATOR_IMPL_ITERATOR_HELPERS_HPP_ */

