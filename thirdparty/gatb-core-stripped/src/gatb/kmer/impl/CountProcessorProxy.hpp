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

#ifndef _COUNT_PROCESSOR_PROXY_HPP_
#define _COUNT_PROCESSOR_PROXY_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/api/ICountProcessor.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

/** */
template<size_t span=KMER_DEFAULT_SPAN>
class CountProcessorProxy : public ICountProcessor<span>
{
public:

    typedef ICountProcessor<span> CountProcessor;
    typedef typename Kmer<span>::Type Type;

    /** Constructor. */
    CountProcessorProxy (ICountProcessor<span>* ref) : _ref(0)  { setRef (ref); }

    /** Destructor. */
    virtual ~CountProcessorProxy()  {  setRef(0); }

    /********************************************************************/
    /*   METHODS CALLED ON THE PROTOTYPE INSTANCE (in the main thread). */
    /********************************************************************/

    /** \copydoc ICountProcessor<span>::begin */
    void begin (const impl::Configuration& config)  { _ref->begin (config); }

    /** \copydoc ICountProcessor<span>::end */
    void end   ()  { _ref->end (); }

    /** \copydoc ICountProcessor<span>::beginPass */
    void beginPass (size_t passId)  { _ref->beginPass (passId); }

    /** \copydoc ICountProcessor<span>::endPass */
    void endPass   (size_t passId)  {  _ref->endPass (passId); }

    /** \copydoc ICountProcessor<span>::clone */
    ICountProcessor<span>* clone ()  { return _ref->clone(); }

    /** \copydoc ICountProcessor<span>::finishClones */
    void finishClones (std::vector<ICountProcessor<span>*>& clones)  { _ref->finishClones (clones); }

    /********************************************************************/
    /*   METHODS CALLED ON ONE CLONED INSTANCE (in a separate thread).  */
    /********************************************************************/

    /** \copydoc ICountProcessor<span>::beginPart */
    void beginPart (size_t passId, size_t partId, size_t cacheSize, const char* name)  { _ref->beginPart (passId, partId, cacheSize, name); }

    /** \copydoc ICountProcessor<span>::endPart */
    void endPart   (size_t passId, size_t partId)  { _ref->endPart (passId, partId); }

    /** \copydoc ICountProcessor<span>::process */
    bool process (size_t partId, const Type& kmer, const CountVector& count, CountNumber sum=0)
    {  return _ref->process (partId, kmer, count, sum);  }

    /*****************************************************************/
    /*                          MISCELLANEOUS.                       */
    /*****************************************************************/

    /** \copydoc ICountProcessor<span>::getName */
    std::string getName() const {  return _ref->getName(); }

    /** \copydoc ICountProcessor<span>::setName */
    void setName (const std::string& name)  { _ref->setName (name); }

    /** \copydoc ICountProcessor<span>::getProperties */
    tools::misc::impl::Properties getProperties() const  { return _ref->getProperties(); }

    /** \copydoc ICountProcessor<span>::getInstances */
    std::vector<ICountProcessor<span>*> getInstances () const  { return _ref->getInstances(); }

protected:

    ICountProcessor<span>* _ref;
    void setRef (ICountProcessor<span>* ref) { SP_SETATTR(ref); }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _COUNT_PROCESSOR_PROXY_HPP_ */
