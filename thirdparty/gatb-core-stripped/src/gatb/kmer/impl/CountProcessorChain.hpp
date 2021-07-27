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

#ifndef _COUNT_PROCESSOR_CHAIN_HPP_
#define _COUNT_PROCESSOR_CHAIN_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/CountProcessorAbstract.hpp>
#include <cstdarg>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

/** The CountProcessorChain implementation allows to link several instances of
 * ICountProcessor. When such an instance is called via 'process', the first item
 * of the list is called via 'process'; if it returns true, the next item in the list
 * is called and so on; if it returns false, the chain is stopped. This class is used
 * for the definition of the "DSK" count processor (histogram -> solidity -> dump)
 */
template<size_t span=KMER_DEFAULT_SPAN>
class CountProcessorChain : public CountProcessorAbstract<span>
{
public:

    typedef ICountProcessor<span> CountProcessor;
    typedef typename Kmer<span>::Type Type;

    /** Constructor. */
    CountProcessorChain (const std::vector<CountProcessor*> items , std::vector<bool>& solidVec) : _items(items) , _solidVec(solidVec) {}

    /** Constructor. */
    CountProcessorChain (CountProcessor* first, ...)
    {
        va_list ap;
        va_start (ap, first);
        for (CountProcessor* loop = first; loop != 0; loop = va_arg(ap, CountProcessor*))
        {
            _items.push_back (loop);  loop->use();
        }
        va_end (ap);
    }

    /** Destructor. */
    virtual ~CountProcessorChain()  {  for (size_t i=0; i<_items.size(); i++)  {  _items[i]->forget(); }  }

    /********************************************************************/
    /*   METHODS CALLED ON THE PROTOTYPE INSTANCE (in the main thread). */
    /********************************************************************/

    /** \copydoc ICountProcessor<span>::begin */
    void begin (const Configuration& config)  {
  		for (size_t i=0; i<_items.size(); i++)  {  _items[i]->begin (config);  }
		this->_solidVec = config._solidVec;

	}

    /** \copydoc ICountProcessor<span>::end */
    void end   ()  {  for (size_t i=0; i<_items.size(); i++)  {  _items[i]->end   ();         }  }

    /** \copydoc ICountProcessor<span>::clone */
    CountProcessorAbstract<span>* clone ()
    {
        std::vector<CountProcessor*> clones;
        for (size_t i=0; i<_items.size(); i++)  {  clones.push_back (_items[i]->clone());  }
        return new CountProcessorChain (clones,_solidVec);
    }

    /** \copydoc ICountProcessor<span>::finishClones */
    void finishClones (std::vector<ICountProcessor<span>*>& clones)
    {
        /** We transform the incoming instances in the correct type. */
        std::vector<CountProcessorChain*> typedClones;
        for (size_t i=0; i<clones.size(); i++)
        {
            if (CountProcessorChain* c = dynamic_cast<CountProcessorChain*>(clones[i]))  { typedClones.push_back(c); }
        }
        if (typedClones.size() != clones.size()) { throw system::Exception("Error in CountProcessorChain::finishClones"); }

        for (size_t i=0; i<_items.size(); i++)
        {
            std::vector<ICountProcessor<span>*> notifyClones (typedClones.size());
            for (size_t j=0; j<typedClones.size(); j++)  {  notifyClones[j] = typedClones[j]->_items[i];  }

            /** We can now notify the N clones for the current ith part in the chain. */
            _items[i]->finishClones (notifyClones);
        }
    }

    /********************************************************************/
    /*   METHODS CALLED ON ONE CLONED INSTANCE (in a separate thread).  */
    /********************************************************************/

    /** \copydoc ICountProcessor<span>::beginPart */
    void beginPart (size_t passId, size_t partId, size_t cacheSize, const char* name)
    {
        for (size_t i=0; i<_items.size(); i++)  {  _items[i]->beginPart (passId, partId, cacheSize, name);  }
    }

    /** \copydoc ICountProcessor<span>::endPart */
    void endPart   (size_t passId, size_t partId)
    {
        for (size_t i=0; i<_items.size(); i++)  {  _items[i]->endPart   (passId, partId);  }
    }

    /** \copydoc ICountProcessor<span>::process */
    bool process (size_t partId, const Type& kmer, const CountVector& count, CountNumber sum=0)
    {
        if (sum==0)  { sum = this->computeSum(count); }
        bool res = true;
        for (size_t i=0; res && i<_items.size(); i++)  {  res = _items[i]->process (partId, kmer, count, sum);  }
        return res;
    }

    /*****************************************************************/
    /*                          MISCELLANEOUS.                       */
    /*****************************************************************/

    /** \copydoc ICountProcessor<span>::getProperties */
    tools::misc::impl::Properties getProperties() const
    {
        tools::misc::impl::Properties result;
        for (size_t i=0; i<_items.size(); i++)  {  result.add (0, _items[i]->getProperties());  }
        return result;
    }

    /** \copydoc ICountProcessor<span>::getInstances */
    std::vector<CountProcessor*> getInstances () const
    {
        std::vector<CountProcessor*> res;
        for (size_t i=0; i<_items.size(); i++)
        {
            std::vector<CountProcessor*> c = _items[i]->getInstances();
            for (size_t i=0; i<c.size(); i++)  { res.push_back(c[i]); }
        }
        return res;
    }

protected:

    /** \copydoc ICountProcessor<span>::computeSum */
    CountNumber computeSum (const CountVector& count) const
    {
        /** Optimization. */
        if (count.size()==1)  { return count[0]; }
        CountNumber sum=0; for (size_t k=0; k<count.size(); k++)  { if (_solidVec.at(k)) sum+=count[k]; }
		return sum;
    }

    std::vector<CountProcessor*> _items;
	
	std::vector<bool> _solidVec;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _COUNT_PROCESSOR_CHAIN_HPP_ */
