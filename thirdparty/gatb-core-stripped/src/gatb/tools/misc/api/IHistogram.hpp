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

/** \file IHistogram.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Interface for histogram (something counting abundances).
 */

#ifndef _GATB_CORE_TOOLS_MISC_IHISTOGRAM_HPP_
#define _GATB_CORE_TOOLS_MISC_IHISTOGRAM_HPP_

#include <gatb/system/api/ISmartPointer.hpp>
#include <gatb/tools/storage/impl/Storage.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
/********************************************************************************/

/** \brief Interface for kmers distribution management
 *
 * This interface allows to have an idea of the function y(x), where x is the occurrence number of a kmer
 * and y is the number of kmers occurring x times.
 *
 * It is often interesting to have a graphical display of this kind of distribution; for instance, it may
 * give an estimation of the coverage of NGS data.
 *
 * We can also find x0 at the first minimum of y(x) : for x<x0, we are likely to have sequencing errors.
 * The first maximum at x1 (x1>x0) is also interesting because it provides an estimation of the reads
 * coverage.
 *
 * This interface is mainly used by the SortingCountAlgorithm.
 *
 * Here is a command line for showing the histogram with gnuplot from the h d f 5 file 'graph.h5'
 *  * h5dump -y -d dsk/histogram graph.h5 | grep [0-9] | grep -v [A-Z].* | paste - - | gnuplot -p -e 'plot [][0:100] "-" with lines'
 *
 *  For the sum of the distribution, you can use:
 *  * h5dump -y -d dsk/histogram graph.h5 | grep [0-9] | grep -v [A-Z].* | paste - - | gawk 'BEGIN{s=0; i=0} { s=s+$2; i=i+1; print i,"  ", s}' | gnuplot -p -e 'plot [0:10][0:] "-" with lines'
*/
class IHistogram : virtual public system::ISmartPointer
{
public:

    /********************************************************************************/
    struct Entry
    {
        u_int16_t index;
        u_int64_t abundance;
        
        /** Comparison operator
         * \param[in] other : object to be compared to
         * \return true if the provided kmer value is greater than the current one. */
        bool operator< (const Entry& other) const {  return this->index < other.index; }
        
        /** Equal operator
         * \param[in] other : object to be compared to
         * \return true if the provided kmer value is greater than the current one. */
        bool operator== (const Entry& other) const {  return (this->index == other.index && this->abundance == other.abundance); }
    };

    /** Destructor. */
    virtual ~IHistogram() {}

    /** Return the maximum allowed for X.
     * \return the max X value. */
    virtual size_t getLength() = 0;

	
	/** Return the maximum allowed for Y in case of 2D histogram.
	 * \return the max Y value. */
	virtual size_t getLength2() = 0;
	
    /** Increase the number of kmers occurring X time
     * \param[in] index : the X value. */
    virtual void inc (u_int16_t index) = 0;

	/** Increase the number of kmers occurring X time in genome and Y times in read
	 * \param[in] index1 : the X value.
	 * \param[in] index2 : the Y value. */
	virtual void inc2D (u_int16_t index1, u_int16_t index2) = 0;
	
	
    /** Save the distribution. It is saved into the bag provided at construction. */
    virtual void save (tools::storage::impl::Group& group) = 0;

	/** Compute first minimum at x0 and firt maximum at x1 (x1>x0). */
    virtual void compute_threshold (int min_auto_threshold) = 0;  //min_auto_threshold = prevents the auto_cutoff from being below this value. Default =3)
	
    /** Get the solid cutoff, ie the x0 at first minimum.
     * \return x0 */
	virtual u_int16_t get_solid_cutoff () = 0;

    /** Get the number of kmers for x>x0, aka solid kmers for x0 threshold
     * \return number of kmers. */
	virtual u_int64_t get_nbsolids_auto () = 0;

	
	/** Get the ratio of weak kmers in total volume
	 * \return ratio */
	virtual float get_ratio_weak () = 0;
	
	/** Get the x1 value at the first maximum after x0. */
    virtual u_int16_t get_first_peak () = 0;

    /** Retrieve the value for x.
     * \param[in] idx : x value.
     * \return y(x). */
    virtual u_int64_t& get (u_int16_t idx) = 0;
	
	
	/** Retrieve the value for x and y of histo2D.
	 * \param[in] idx1 : x value.
	 * \param[in] idx2 : y value.
	 * \return cpt(x,y). */
	virtual u_int64_t& get2D (u_int16_t idx1,u_int16_t idx2) = 0;
	
	
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_IHISTOGRAM_HPP_ */
