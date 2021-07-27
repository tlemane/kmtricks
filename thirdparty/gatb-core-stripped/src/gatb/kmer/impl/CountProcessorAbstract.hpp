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

#ifndef _COUNT_PROCESSOR_ABSTRACT_HPP_
#define _COUNT_PROCESSOR_ABSTRACT_HPP_

/********************************************************************************/

#include <gatb/kmer/api/ICountProcessor.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

/** Abstract implementation of ICountProcessor interface.
 */
template<size_t span>
class CountProcessorAbstract : public ICountProcessor<span>
{
public:

    /** Shortcuts. */
    typedef typename Kmer<span>::Type Type;

    /** Constructor. */
    CountProcessorAbstract (const std::string& name="processor") : _name(name)  {}

    /** Destructor. */
    virtual ~CountProcessorAbstract ()  {}

    /********************************************************************/
    /*   METHODS CALLED ON THE PROTOTYPE INSTANCE (in the main thread). */
    /********************************************************************/

    /** \copydoc ICountProcessor<span>::begin */
    virtual void begin(const Configuration& config) {}

    /** \copydoc ICountProcessor<span>::end */
    virtual void end  () {}

    /** \copydoc ICountProcessor<span>::end */
    virtual void beginPass (size_t passId) {}

    /** \copydoc ICountProcessor<span>::end */
    virtual void endPass   (size_t passId) {}

    /** \copydoc ICountProcessor<span>::finishClones */
    virtual void finishClones (std::vector<ICountProcessor<span>*>& clones)  {}

    /********************************************************************/
    /*   METHODS CALLED ON ONE CLONED INSTANCE (in a separate thread).  */
    /********************************************************************/

    /** \copydoc ICountProcessor<span>::beginPart */
    virtual void beginPart (size_t passId, size_t partId, size_t cacheSize, const char* name) {}

    /** \copydoc ICountProcessor<span>::endPart */
    virtual void endPart   (size_t passId, size_t partId) {}

    /** \copydoc ICountProcessor<span>::process */
    virtual bool process (size_t partId, const Type& kmer, const CountVector& count, CountNumber sum=0)  {  return true;  }

    /*****************************************************************/
    /*                          MISCELLANEOUS.                       */
    /*****************************************************************/

    /** \copydoc ICountProcessor<span>::getName */
    virtual std::string getName() const  { return _name; }

    /** \copydoc ICountProcessor<span>::setName */
    virtual void setName (const std::string& name) { _name = name; }

    /** \copydoc ICountProcessor<span>::getProperties */
    virtual tools::misc::impl::Properties getProperties() const { return tools::misc::impl::Properties(); }

    /** \copydoc ICountProcessor<span>::getInstances */
    virtual std::vector<ICountProcessor<span>*> getInstances () const
    {
        std::vector<ICountProcessor<span>*> res;
        res.push_back((CountProcessorAbstract*)this);
        return res;
    }

private:

    std::string _name;
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _COUNT_PROCESSOR_ABSTRACT_HPP_ */
