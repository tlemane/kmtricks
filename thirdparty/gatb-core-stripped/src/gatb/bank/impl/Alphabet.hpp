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

/** \file Alphabet.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Implementation for genomic alphabets
 */

#ifndef _GATB_CORE_BANK__IMPL_ALPHABET_HPP_
#define _GATB_CORE_BANK__IMPL_ALPHABET_HPP_

/********************************************************************************/

#include <gatb/bank/api/IAlphabet.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for genomic databases management. */
namespace bank      {
/** \brief Implementation for genomic databases management. */
namespace impl      {
/********************************************************************************/

/* Define the kind of the underlying alphabet. For instance, it could be the kind of sequences
 * (protein, ADN...) read from a FASTA file
 */
class AlphabetNucleic : public IAlphabet
{
public:

    /** Singleton method.
     * \return the singleton instance. */
    static IAlphabet& singleton()  { static AlphabetNucleic instance; return instance; }

    /** \copydoc IAlphabet::getKind */
    Kind_e getKind ()  { return NUCLEIC_ACID; }

    /** \copydoc IAlphabet::getLetters */
    std::string  getLetters ()  { return std::string ("ACTG"); }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_BANK__IMPL_ALPHABET_HPP_ */
