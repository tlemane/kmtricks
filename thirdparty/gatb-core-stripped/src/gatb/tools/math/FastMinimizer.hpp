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

/** \file FastMinimizer.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Fast computation of lexicographical minimizers wtih no-AA-inside constraint
 */


#ifndef _GATB_CORE_TOOLS_MATH_FASTMINIMIZER_HPP_
#define _GATB_CORE_TOOLS_MATH_FASTMINIMIZER_HPP_


#include <stdint.h>

extern const unsigned char revcomp_4NT[];

template<typename T, typename minimizer_type> inline void fastLexiMinimizerChunk (T val, const unsigned int _nbMinimizers, const unsigned int m, const minimizer_type high_bits, minimizer_type &minimizer, size_t &position, size_t position_offset, bool &AA_found) 
{
    // FIXME: useless for minimizer size larger than 16 just because of the 0x55555555 mask 
    if (m > 16) {AA_found = false; return;}
    
    /* those require only a single AND operation, rest are constants */
    #define BINARY_MMER_STARTS_WITH_2MER(val,m,binnucl) ((val & (15 << (2*(m-2)))) == (binnucl << (2*(m-2))) )
    #define BINARY_MMER_ENDS_WITH_2MER(val,binnucl)  ((val & 15) == binnucl)
    unsigned int j = 0;
    bool adjusted_near_end = (high_bits == 0); /* will adjust only if there are bits */
    const unsigned int nb_minim_in_chunk = sizeof(T)*4;
    const minimizer_type m_mask = (1 << (2*m)) - 1;
    
    //int next_TT_in_shifts = -1;
 
    /* numbers of minimizers to examine in this chunk
     * we assume that each nucleotide of a chunk is the start of a minimizer
     * because we add missing (m-1) nucl as high bits later */
    const unsigned int it = std::min((unsigned int)nb_minim_in_chunk, _nbMinimizers - (unsigned int)position_offset); 

    while (j < it)
    {
        if ((j >= nb_minim_in_chunk - m) && (! adjusted_near_end))
        { 
            /* append the next m-1 (or less) characters to "val", in order 
             * to keep iterating minimizers smoothly across the boundary of our representation*/
            val |= high_bits << ((nb_minim_in_chunk - j)*2);
            adjusted_near_end = true;
        }

        //bool mmer_starts_with_TT = (BINARY_MMER_STARTS_WITH_2MER(val,m,10));
        bool mmer_starts_with_AA = (BINARY_MMER_STARTS_WITH_2MER(val,m,0  /* (AA=0b) */));
        bool mmer_ends_with_TT = (BINARY_MMER_ENDS_WITH_2MER(val,10 /* (TT = 1010b = 10)*/));

        minimizer_type candidate = val & m_mask; 

        if (mmer_ends_with_TT)
        {
            minimizer_type candidate_revcomp = ((revcomp_4NT [val&0xFF] << 8) | revcomp_4NT [(val>>8)&0xFF]) & m_mask; 
            if (mmer_starts_with_AA) 
                candidate = std::min(candidate, candidate_revcomp);
            else
                candidate = candidate_revcomp;
        }

        if (mmer_starts_with_AA || mmer_ends_with_TT)
        {
            AA_found = true;

            if (candidate < minimizer)
            {
                /* check if it's allowed
                   for comments of that code, see the is_allowed function */
                minimizer_type a1 = candidate & m_mask; //
                a1 =   ~(( a1 )   | (  a1 >>2 ));  //
                a1 =((a1 >>1) & a1) & ((1 << ((m-2)*2)) -1 ) & 0x55555555;  //
                
                if (a1 == 0)
                {
                    /* set new minimizer */
                    minimizer = candidate;

                    /* last minimizer has position k - m 
                     * nbminimizers is k - m 
                     * position_offset is here because k is divided by chunks of 32 nucleotides */
                    position = (_nbMinimizers - 1) - position_offset - j;
                }
            }

        } // end if (TT last or AA first)
        // commented, because turns out this was slower in practice (k=51, e.coli).
        /*

         * let's call AA-type a m-mer of the form AA......
         * and TT-type a m-mer of the form ......TT
         * (it is the reverse-complement of a possible minimizer of the kmer)
         * m-mers which are both AA-type and TT-type are corner cases.
         *
         * furthermore, we only care about AA-type or TT-type m-mers ("candidates"). 
         * other m-mers are ignored.
         *
         * we make the hypothesis (see above) that the minimizer is either a candidate 
         * (if it's AA-type) or the reverse complement of a candidate (if it's TT-type).
         *
         * idea: a minimizer can never be of the form AA...AA....
         * or ....TT....TT, so, when we see AA......, we know that we can right shift
         * (since we read the binary kmer from end to beginning) m-1 positions before *
         * getting another AA-type candidate.
         * however, doing so, we might miss a TT-type candidate.
         *
         * so let's keep track of the TT...... m-mers. whenever we haven't seen one for a while, 
         * means that we can skip safely.

         if (mmer_starts_with_TT)
         {
             if (next_TT_in_shifts == -1)
                 next_TT_in_shifts = m-2;
             else
                 next_TT_in_shifts = std::min(next_TT_in_shifts,m-2);
         }

         if (mmer_starts_with_AA && next_TT_in_shifts != -1 && j >= m)
         {
             val >>= 12;
             j += 6;
             next_TT_in_shifts = -1;
             continue;
         }

         next_TT_in_shifts = std::max(next_TT_in_shifts-1,-1);
         */


        val >>= 2;
        j++;

    } // end while (minimizers)
}

#endif
