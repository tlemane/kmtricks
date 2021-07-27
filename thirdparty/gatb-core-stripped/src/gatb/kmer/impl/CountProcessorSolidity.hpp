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

#ifndef _COUNT_PROCESSOR_SOLIDITY_HPP_
#define _COUNT_PROCESSOR_SOLIDITY_HPP_

/********************************************************************************/

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/CountProcessorAbstract.hpp>
#include <gatb/tools/misc/api/Enums.hpp>
#include <gatb/tools/misc/api/Range.hpp>
#include <gatb/tools/misc/impl/Stringify.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

class CountProcessorSolidityInfo
{
public:
    CountProcessorSolidityInfo () {}
	CountProcessorSolidityInfo (const std::vector<tools::misc::CountRange>& thresholds, std::vector<bool> solidVec) : _thresholds(thresholds), _solidVec(solidVec) {};

    /** Update abundance min in the threshold ranges. */
    void setAbundanceMin (const std::vector<CountNumber>& cutoffs)
    {
        if (cutoffs.size() > _thresholds.size())
        {
            throw system::Exception ("Unable to set abundance min values (%d values for %d banks)", cutoffs.size(), _thresholds.size());
        }

        std::vector<tools::misc::CountRange> newThresholds;
        for (size_t i=0; i<cutoffs.size(); i++)
        {
            /** We change the value only if it is "auto". */
            CountNumber abundanceMin = _thresholds[i].getBegin()==-1 ? cutoffs[i] :  _thresholds[i].getBegin();
            newThresholds.push_back (tools::misc::CountRange(abundanceMin, _thresholds[i].getEnd() ) );
        }

        /** We may have less cutoffs than banks; we complete with the last cutoff value. */
        tools::misc::CountRange lastRange (newThresholds[newThresholds.size()-1]);
        for (size_t i=cutoffs.size(); i<_thresholds.size(); i++)  {  newThresholds.push_back (lastRange);  }

        /** We update the thresholds for the current count processor. */
        _thresholds = newThresholds;
    }

protected:
    std::vector<tools::misc::CountRange> _thresholds;
	std::vector<bool> _solidVec;
};

/********************************************************************************/

/** The CountProcessorSolidityAbstract is an abstract class that factories stuff
 * for telling whether a kmer is solid or not.
 *
 * Inherited classes provides (through the 'check' method) the way the kmer solidity is
 * computed. There is one subclass per kind of kmer solidity.
 *
 * Technically, static polymorphism is used here through the 'check' method.
 *
 * Note that there exists a factory class CountProcessorSolidityFactory that manages
 * the creation of the correct instance according to some user information.
 */
template<size_t span, class Derived>
class CountProcessorSolidityAbstract : public CountProcessorAbstract<span>, public CountProcessorSolidityInfo
{
public:

    /** Constructor for prototype instance. */
    CountProcessorSolidityAbstract () : _total(0), _ok(0)   {}

    /** Constructor for clone instance. */
	CountProcessorSolidityAbstract (const std::vector<tools::misc::CountRange>& thresholds, std::vector<bool>& solidVec)
        : CountProcessorSolidityInfo(thresholds,solidVec), _total(0), _ok(0)   {}

    /** Destructor. */
    virtual ~CountProcessorSolidityAbstract()  {}

    /********************************************************************/
    /*   METHODS CALLED ON THE PROTOTYPE INSTANCE (in the main thread). */
    /********************************************************************/

    /** \copydoc ICountProcessor<span>::begin */
    void begin(const Configuration& config)
    {
        /** We copy the abundance thresholds got during configuration. */
        this->_thresholds = config._abundance;
		this->_solidVec = config._solidVec;
    }

    /** \copydoc ICountProcessor<span>::clones */
    CountProcessorAbstract<span>* clone ()  { return new Derived (_thresholds,_solidVec); }

    /** \copydoc ICountProcessor<span>::finishClones */
    void finishClones (std::vector<ICountProcessor<span>*>& clones)
    {
        for (size_t i=0; i<clones.size(); i++)
        {
            /** We have to recover type information. */
            if (CountProcessorSolidityAbstract* clone = dynamic_cast<CountProcessorSolidityAbstract*> (clones[i]))
            {
                this->_total += clone->_total;
                this->_ok    += clone->_ok;
            }
        }
    }

    /********************************************************************/
    /*   METHODS CALLED ON ONE CLONED INSTANCE (in a separate thread).  */
    /********************************************************************/

    /** \copydoc ICountProcessor<span>::process */
    bool process (size_t partId, const typename Kmer<span>::Type& kmer, const CountVector& count, CountNumber sum)
    {
        /** We use static polymorphism here. */
        bool result = static_cast<Derived*>(this)->check (count, sum);

        _total ++;
        if (result)  { _ok++; }
        return result;
    }

    /*****************************************************************/
    /*                          MISCELLANEOUS.                       */
    /*****************************************************************/

    /** \copydoc ICountProcessor<span>::getProperties */
    tools::misc::impl::Properties getProperties() const
    {
        tools::misc::impl::Properties result;
        result.add (0, "kmers");
        result.add (1, "solidity_kind",      "%s", this->getName().c_str());

        std::stringstream ss;
        for (size_t i=0; i<_thresholds.size(); i++)  {  ss << _thresholds[i].getBegin() << " "; }
        result.add (1, "thresholds",         ss.str().c_str());

        result.add (1, "kmers_nb_distinct",  "%ld", _total);
        result.add (1, "kmers_nb_solid",     "%ld", _ok);
        result.add (1, "kmers_nb_weak",      "%ld", _total - _ok);
        if (_total > 0)  {  result.add (1, "kmers_percent_weak", "%.1f", 100.0 - 100.0 * (double)_ok / (double)_total  );  }

        return result;
    }

protected:

    u_int64_t _total;
    u_int64_t _ok;
};

/********************************************************************************/

template<size_t span>
class CountProcessorSoliditySum : public CountProcessorSolidityAbstract<span,CountProcessorSoliditySum<span> >
{
public:

    CountProcessorSoliditySum () {}

    CountProcessorSoliditySum (const std::vector<tools::misc::CountRange>& thresholds, std::vector<bool>& solidVec)
        : CountProcessorSolidityAbstract<span,CountProcessorSoliditySum<span> > (thresholds,solidVec)  {}

    bool check (const CountVector& count, CountNumber sum)
    {
        return this->_thresholds[0].includes (sum);
    }

    std::string getName() const  { return std::string("sum"); }
};

/********************************************************************************/

template<size_t span>
class CountProcessorSolidityMax : public CountProcessorSolidityAbstract<span,CountProcessorSolidityMax<span> >
{
public:

    CountProcessorSolidityMax () {}

	CountProcessorSolidityMax (const std::vector<tools::misc::CountRange>& thresholds, std::vector<bool>& solidVec)
        : CountProcessorSolidityAbstract<span,CountProcessorSolidityMax<span> > (thresholds,solidVec)  {}

    bool check (const CountVector& count, CountNumber sum)
    {
        return this->_thresholds[0].includes (*std::max_element (count.begin(),count.end()));
    }

    std::string getName() const  { return std::string("max"); }
};

/********************************************************************************/

template<size_t span>
class CountProcessorSolidityMin : public CountProcessorSolidityAbstract<span,CountProcessorSolidityMin<span> >
{
public:

    CountProcessorSolidityMin () {}

    CountProcessorSolidityMin (const std::vector<tools::misc::CountRange>& thresholds, std::vector<bool>& solidVec)
        : CountProcessorSolidityAbstract<span,CountProcessorSolidityMin<span> > (thresholds,solidVec)  {}

    bool check (const CountVector& count, CountNumber sum)
    {
        return this->_thresholds[0].includes (*std::min_element (count.begin(),count.end()));
    }

    std::string getName() const  { return std::string("min"); }
};

/********************************************************************************/

template<size_t span>
class CountProcessorSolidityAll : public CountProcessorSolidityAbstract<span,CountProcessorSolidityAll<span> >
{
public:

    CountProcessorSolidityAll () {}

    CountProcessorSolidityAll (const std::vector<tools::misc::CountRange>& thresholds, std::vector<bool>& solidVec)
        : CountProcessorSolidityAbstract<span,CountProcessorSolidityAll<span> > (thresholds,solidVec)  {}

    bool check (const CountVector& count, CountNumber sum)
    {
        for (size_t i=0; i<count.size(); i++)  {  if (this->_thresholds[i].includes(count[i]) == false)   { return false; }  }
        return true;
    }

    std::string getName() const  { return std::string("all"); }
};

/********************************************************************************/

template<size_t span>
class CountProcessorSolidityOne : public CountProcessorSolidityAbstract<span, CountProcessorSolidityOne<span> >
{
public:

    CountProcessorSolidityOne () {}

    CountProcessorSolidityOne (const std::vector<tools::misc::CountRange>& thresholds, std::vector<bool>& solidVec)
        : CountProcessorSolidityAbstract<span, CountProcessorSolidityOne<span> > (thresholds,solidVec)  {}

    bool check (const CountVector& count, CountNumber sum)
    {
        for (size_t i=0; i<count.size(); i++)  {  if (this->_thresholds[i].includes(count[i]) == true)   { return true; }  }
        return false;
    }

    std::string getName() const  { return std::string("one"); }
};
	

/********************************************************************************/

	
	template<size_t span>
	class CountProcessorSolidityCustom : public CountProcessorSolidityAbstract<span, CountProcessorSolidityCustom<span> >
	{
	public:
		
		CountProcessorSolidityCustom () {}
		
		CountProcessorSolidityCustom (const std::vector<tools::misc::CountRange>& thresholds, std::vector<bool>& solidVec)
		: CountProcessorSolidityAbstract<span, CountProcessorSolidityCustom<span> > (thresholds,solidVec)  {}
		
		
		bool check (const CountVector& count, CountNumber sum)
		{
			for (size_t i=0; i<count.size(); i++)  {
				
				if (this->_solidVec.at(i) == false &&   this->_thresholds[i].includes(count[i]) == true   )   { return false; }
				else if (this->_solidVec.at(i) == true &&   this->_thresholds[i].includes(count[i]) == false  ) { return false; }
			
			}
			return true;
		}
		
		std::string getName() const  { return std::string("custom"); }
		
	};

	

/********************************************************************************/

template<size_t span>
class CountProcessorSolidityFactory
{
public:

    /** */
	static ICountProcessor<span>* create (tools::misc::KmerSolidityKind kind)
    {
        switch (kind)
        {
        case tools::misc::KMER_SOLIDITY_MIN: return new CountProcessorSolidityMin<span> ();
        case tools::misc::KMER_SOLIDITY_MAX: return new CountProcessorSolidityMax<span> ();
        case tools::misc::KMER_SOLIDITY_ONE: return new CountProcessorSolidityOne<span> ();
		case tools::misc::KMER_SOLIDITY_CUSTOM: return new CountProcessorSolidityCustom<span> ();
        case tools::misc::KMER_SOLIDITY_ALL: return new CountProcessorSolidityAll<span> ();
        case tools::misc::KMER_SOLIDITY_SUM: return new CountProcessorSoliditySum<span> ();
        default:  throw system::Exception ("unable to create CountProcessorSolidity instance for kind %d", kind);
        }
    }

    /** */
    static ICountProcessor<span>* create (tools::misc::IProperties& props)
    {
        tools::misc::KmerSolidityKind kind;
        parse (props.getStr (STR_SOLIDITY_KIND), kind);
		return create (kind);
    }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _COUNT_PROCESSOR_SOLIDITY_HPP_ */
