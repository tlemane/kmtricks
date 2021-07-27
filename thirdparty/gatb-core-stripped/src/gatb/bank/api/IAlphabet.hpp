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

/** \file IAlphabet.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Interface definition for genomic alphabets
 */

#ifndef _GATB_CORE_BANK_IALPHABET_HPP_
#define _GATB_CORE_BANK_IALPHABET_HPP_

/********************************************************************************/

#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for genomic databases management. */
namespace bank      {
/********************************************************************************/

/* Define the kind of the underlying alphabet. For instance, it could be the kind of sequences
 * (protein, ADN...) read from a FASTA file.
 */
class IAlphabet
{
public:

    /** Define the different kinds of alphabets. */
    enum Kind_e
    {
        /** Amino acid alphabet */
        AMINO_ACID,
        /** Nucleic acid alphabet */
        NUCLEIC_ACID
    };

    /** \return the kind of the alphabet. */
    virtual Kind_e getKind () = 0;

    /** Return the letters that make the alphabet. For instance, the nucleic alphabet should return "ACGT".
     * \return the alphabet letters.*/
    virtual std::string  getLetters () = 0;

    /** Destructor. */
    virtual ~IAlphabet () {}
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK_IALPHABET_HPP_ */
