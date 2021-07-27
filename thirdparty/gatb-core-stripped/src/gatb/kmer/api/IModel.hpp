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

/** \file IModel.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Interface definition for the kmer model
 */

#ifndef _GATB_CORE_KMER_IMODEL_HPP_
#define _GATB_CORE_KMER_IMODEL_HPP_

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
/********************************************************************************/

/** \brief enumeration giving the way the kmers are computed. */
enum KmerMode
{
    /** Kmer from direct strand */
    KMER_DIRECT,
    /** Kmer as reverse complement */
    KMER_REVCOMP,
    /** Kmer as minimum between the direct and reverse complement */
    KMER_CANONICAL
};

/********************************************************************************/

/** Enumeration for strand. */
enum Strand
{
    /** enum for the forward strand */
    STRAND_FORWARD = (1<<0),
    /** enum for the reverse strand */
    STRAND_REVCOMP = (1<<1),
    STRAND_ALL     = STRAND_FORWARD + STRAND_REVCOMP
};

/** Ascii representation of a Strand object  ("FWD" or "REV")
 * \param[in] s : the strand
 * \return the Ascii representation */
inline std::string toString (Strand s)
{
         if (s == STRAND_FORWARD)  { return std::string("FWD"); }
    else if (s == STRAND_REVCOMP)  { return std::string("REV"); }
    else { return std::string("???"); }
}

inline Strand StrandReverse (const Strand& s)  {  return (s==STRAND_FORWARD ? STRAND_REVCOMP : STRAND_FORWARD);  }

/********************************************************************************/

/** Enumeration for nucleotides. */
enum Nucleotide
{
    /** enum for the A nucleotide */
    NUCL_A   = 0,
    /** enum for the C nucleotide */
    NUCL_C   = 1,
    /** enum for the T nucleotide */
    NUCL_T   = 2,
    /** enum for the G nucleotide */
    NUCL_G   = 3,
    NUCL_UNKNOWN = 4
};

/** Ascii representation of a Nucleotide object ('A', 'C', 'T' or 'G')
 * \param[in] nt : the nucleotide
 * \return the Ascii representation */
inline char ascii (Nucleotide nt)
{
    static char table[] = {'A', 'C', 'T', 'G', 'N' };
    return table[(int)nt];
}

/** Reverse complement of a nucleotide (A->T, C->G, T->A, G->C)
 * \param[in] nt : the nucleotide to be reversed
 * \return the reversed nucleotide
 */
inline Nucleotide reverse (Nucleotide nt)
{
    static Nucleotide table[] = {NUCL_T, NUCL_G, NUCL_A, NUCL_C, NUCL_UNKNOWN};
    return table[(int)nt];
}

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMODEL_HPP_ */
