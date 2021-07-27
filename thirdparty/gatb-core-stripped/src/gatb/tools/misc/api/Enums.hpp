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

/** \file Enums.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Enumerations
 */

#ifndef _GATB_CORE_TOOLS_MISC_ENUMS_HPP_
#define _GATB_CORE_TOOLS_MISC_ENUMS_HPP_

#include <gatb/system/api/Exception.hpp>
#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
/********************************************************************************/

/** Enumeration for the different kinds of bank conversions supported in GATB. */
enum BankConvertKind  {  BANK_CONVERT_NONE, BANK_CONVERT_TMP, BANK_CONVERT_KEEP };

/** Get the enum from a string.
 * \param[in] s : string to be parsed
 * \param[out] kind : enum to be set from the string parsing. */
static void parse (const std::string& s, BankConvertKind& kind)
{
         if (s == "none")   { kind = BANK_CONVERT_NONE;  }
    else if (s == "tmp")    { kind = BANK_CONVERT_TMP;  }
    else if (s == "keep")   { kind = BANK_CONVERT_KEEP;  }
    else   { throw system::Exception ("bad bank convert kind '%s'", s.c_str()); }
}

/** Get the string associated to an enum
 * \param[in] kind : the enum value
 * \return the associated string */
static const char* toString (BankConvertKind kind)
{
    switch (kind)
    {
        case BANK_CONVERT_NONE:     return "none";
        case BANK_CONVERT_TMP:      return "tmp";
        case BANK_CONVERT_KEEP:     return "keep";
        default:        throw system::Exception ("bad bank convert kind %d", kind);
    }
}

/********************************************************************************/

/** Enumeration for the different kinds of Bloom filters supported in GATB. */
enum BloomKind
{
    /** No Bloom filter */
    BLOOM_NONE,
    /** Trivial implementation of Bloom filters */
    BLOOM_BASIC,
    /** Implementation of Bloom filters improving CPU cache management. */
    BLOOM_CACHE,
    /** Implementation of Bloom filters improving CPU cache management. */
    BLOOM_NEIGHBOR,
    BLOOM_DEFAULT
};

/** Get the enum from a string.
 * \param[in] s : string to be parsed
 * \param[out] kind : enum to be set from the string parsing. */
static void parse (const std::string& s, BloomKind& kind)
{
         if (s == "none")        { kind = BLOOM_NONE;  }
    else if (s == "basic")       { kind = BLOOM_BASIC;  }
    else if (s == "cache")       { kind = BLOOM_CACHE; }
	else if (s == "neighbor")    { kind = BLOOM_NEIGHBOR; }
    else if (s == "default")     { kind = BLOOM_CACHE; }
    else   { throw system::Exception ("bad Bloom kind '%s'", s.c_str()); }
}

/** Get the string associated to an enum
 * \param[in] kind : the enum value
 * \return the associated string */
static const char* toString (BloomKind kind)
{
    switch (kind)
    {
        case BLOOM_NONE:      return "none";
        case BLOOM_BASIC:     return "basic";
        case BLOOM_CACHE:     return "cache";
		case BLOOM_NEIGHBOR:  return "neighbor";
        case BLOOM_DEFAULT:   return "cache";
        default:        throw system::Exception ("bad Bloom kind %d", kind);
    }
}

/********************************************************************************/

/** Enumeration for the different kinds of cFP storage mechanisms supported in GATB. */
enum DebloomKind
{
    /** No cFP */
    DEBLOOM_NONE,
    /** Save cFP in the original way (a sorted vector) */
    DEBLOOM_ORIGINAL,
    /** Save cFP with cascading Bloom filters. */
    DEBLOOM_CASCADING,
    DEBLOOM_DEFAULT
};

/** Get the debloom kind from a string.
 * \param[in] s : string to be parsed
 * \param[out] kind : the debloom kind to be set from the string parsing. */
static void parse (const std::string& s, DebloomKind& kind)
{
         if (s == "none")       { kind = DEBLOOM_NONE;      }
    else if (s == "original")   { kind = DEBLOOM_ORIGINAL;  }
    else if (s == "cascading")  { kind = DEBLOOM_CASCADING; }
    else if (s == "default")    { kind = DEBLOOM_CASCADING; }
    else   { throw system::Exception ("bad debloom kind '%s'", s.c_str()); }
}

/** Get the string associated to an enum
 * \param[in] kind : the enum value
 * \return the associated string */
static std::string toString (DebloomKind kind)
{
    switch (kind)
    {
        case DEBLOOM_NONE:      return "none";
        case DEBLOOM_ORIGINAL:  return "original";
        case DEBLOOM_CASCADING: return "cascading";
        case DEBLOOM_DEFAULT:   return "cascading";
        default:        throw system::Exception ("bad debloom kind %d", kind);
    }
}

/********************************************************************************/

/** Enumeration for the different kinds of debloom algorithms supported in GATB. */
enum DebloomImpl
{
    /** Initial debloom algorithm */
    DEBLOOM_IMPL_BASIC,
    /** Debloom algorithm based on minimizers. */
    DEBLOOM_IMPL_MINIMIZER,
    DEBLOOM_IMPL_DEFAULT
};

/** Get the enum from a string.
 * \param[in] s : string to be parsed
 * \param[out] kind : enum to be set from the string parsing. */
static void parse (const std::string& s, DebloomImpl& kind)
{
         if (s == "basic")       { kind = DEBLOOM_IMPL_BASIC;      }
    else if (s == "minimizer")   { kind = DEBLOOM_IMPL_MINIMIZER;  }
    else if (s == "default")     { kind = DEBLOOM_IMPL_MINIMIZER; }
    else   { throw system::Exception ("bad debloom impl '%s'", s.c_str()); }
}

/** Get the string associated to an enum
 * \param[in] kind : the enum value
 * \return the associated string */
static std::string toString (DebloomImpl kind)
{
    switch (kind)
    {
        case DEBLOOM_IMPL_BASIC:      return "basic";
        case DEBLOOM_IMPL_MINIMIZER:  return "minimizer";
        case DEBLOOM_IMPL_DEFAULT:    return "minimizer";
        default:        throw system::Exception ("bad debloom impl %d", kind);
    }
}

/********************************************************************************/

/** Enumeration for the different kinds of branching storages supported in GATB. */
enum BranchingKind
{
    /** Don't save branching nodes. */
    BRANCHING_NONE,
    /** Save branching nodes within the graph file. */
    BRANCHING_STORED
};

/** Get the enum from a string.
 * \param[in] s : string to be parsed
 * \param[out] kind : enum to be set from the string parsing. */
static void parse (const std::string& s, BranchingKind& kind)
{
         if (s == "none")     { kind = BRANCHING_NONE;  }
    else if (s == "stored")   { kind = BRANCHING_STORED;  }
    else   { throw system::Exception ("bad branching kind '%s'", s.c_str()); }
}

/** Get the string associated to an enum
 * \param[in] kind : the enum value
 * \return the associated string */
static std::string toString (BranchingKind kind)
{
    switch (kind)
    {
        case BRANCHING_NONE:     return "none";
        case BRANCHING_STORED:   return "stored";
        default:        throw system::Exception ("bad branching kind %d", kind);
    }
}

/********************************************************************************/

/** Enumeration for the different kinds of kmer solidity criteria supported in GATB. */
enum KmerSolidityKind
{
    /** min criteria */
    KMER_SOLIDITY_MIN,
    /** max criteria */
    KMER_SOLIDITY_MAX,
    /** on criteria */
    KMER_SOLIDITY_ONE,
	/** custom criteria */
	KMER_SOLIDITY_CUSTOM,
    /** all criteria */
    KMER_SOLIDITY_ALL,
    /** sum criteria */
    KMER_SOLIDITY_SUM,
    KMER_SOLIDITY_DEFAULT
};

/** Get the enum from a string.
 * \param[in] s : string to be parsed
 * \param[out] kind : enum to be set from the string parsing. */
static void parse (const std::string& s, KmerSolidityKind& kind)
{
         if (s == "min")    { kind = KMER_SOLIDITY_MIN;  }
    else if (s == "max")    { kind = KMER_SOLIDITY_MAX;  }
    else if (s == "one")    { kind = KMER_SOLIDITY_ONE;  }
    else if (s == "all")    { kind = KMER_SOLIDITY_ALL;  }
    else if (s == "sum")    { kind = KMER_SOLIDITY_SUM;  }
	else if ( s.find("custom") != std::string::npos ) { kind = KMER_SOLIDITY_CUSTOM; }
    else   { throw system::Exception ("bad kmer solidity kind '%s'", s.c_str()); }
}

/** Get the string associated to an enum
 * \param[in] kind : the enum value
 * \return the associated string */
static std::string toString (KmerSolidityKind kind)
{
    switch (kind)
    {
        case KMER_SOLIDITY_MIN:     return "min";
        case KMER_SOLIDITY_MAX:     return "max";
        case KMER_SOLIDITY_ONE:     return "one";
		case	KMER_SOLIDITY_CUSTOM: return "custom";
        case KMER_SOLIDITY_ALL:     return "all";
        case KMER_SOLIDITY_SUM:     return "sum";
        case KMER_SOLIDITY_DEFAULT: return "sum";
        default:    throw system::Exception ("bad kmer solidity kind %d", kind);
    }
}

/********************************************************************************/

/** Enumeration of different kinds of graph traversal. */
enum TraversalKind
{
    /** Undefined */
    TRAVERSAL_NONE=0,
    /** Path are unitigs */
    TRAVERSAL_UNITIG=1,
    /** Path are contigs */
    TRAVERSAL_CONTIG=2
};

/** Get the enum from a string.
 * \param[in] s : string to be parsed
 * \param[out] kind : enum to be set from the string parsing. */
static void parse (const std::string& s, TraversalKind& kind)
{
         if (s == "none")      { kind = TRAVERSAL_NONE;  }
    else if (s == "unitig")    { kind = TRAVERSAL_UNITIG;  }
    else if (s == "contig")    { kind = TRAVERSAL_CONTIG;  }
    else   { throw system::Exception ("bad traversal kind '%s'", s.c_str()); }
}

/** Get the string associated to an enum
 * \param[in] kind : the enum value
 * \return the associated string */
static std::string toString (TraversalKind kind)
{
    switch (kind)
    {
        case TRAVERSAL_NONE:    return "none";
        case TRAVERSAL_UNITIG:  return "unitig";
        case TRAVERSAL_CONTIG:  return "contig";
        default:    throw system::Exception ("bad traversal kind %d", kind);
    }
}

/********************************************************************************/

/** Provide different modes for graph traversal stop criteria. */
enum ExtendStopMode_e
{
    /** Stop traversal after the first unitig/contig. */
    ExtendStopMode_after_first_contig,
    /** Stop traversal when the maximum depth is reached. */
    ExtendStopMode_until_max_depth
};

/********************************************************************************/

/** Provide different modes of graph traversal for building extensions. */
enum SearchMode_e
{
    /** Breadth first traversal. */
    SearchMode_Breadth,
    /** Depth first traversal. */
    SearchMode_Depth
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_ENUMS_HPP_ */
