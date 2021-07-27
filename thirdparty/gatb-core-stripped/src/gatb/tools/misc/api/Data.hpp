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

/** \file Data.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Data structure
 */

#ifndef _GATB_CORE_TOOLS_MISC_DATA_HPP_
#define _GATB_CORE_TOOLS_MISC_DATA_HPP_

/********************************************************************************/

#include <gatb/system/api/ISmartPointer.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/api/Range.hpp>
#include <gatb/tools/misc/api/Vector.hpp>

#include <iostream>
#include <string.h>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
/********************************************************************************/

/** \brief Definition of a data chunk
 *
 * A data is defined by:
 *      - an encoding format
 *      - a buffer holding the actual data
 *      - the size of the data
 *
 * It is implemented as a subclass of the Vector class, which allows to define Data
 * as a sub part of a referred Data instance.
 *
 * For instance, Data is used for storing nucleotides sequences inside the Sequence
 * structure.
 *
 * \note
 * In contrast to Vector, the size() represents the number of data elements,
 * not the number of bytes (see Data::getBufferLength() for length in bytes).
 */
class Data : public Vector<char>
{
public:

    /** Define how data is encoded. */
    enum Encoding_e
    {
        /** data encoded as ASCII codes (so one byte per data unit) */
        ASCII,
        /** one byte per data as integer value (for instance: A=0, C=1, T=2, G=3) */
        INTEGER,
        /** 4 nucleotides compressed in one byte */
        BINARY
    };

    /** Default constructor. */
    Data (Encoding_e encode = BINARY)  : encoding(encode) {}

    /** Default constructor. */
    Data (char* buffer)  : encoding(ASCII) { setRef(buffer,strlen(buffer)); }

    /** Constructor. */
    Data (size_t len, Encoding_e encode = BINARY)  : Vector<char>(len), encoding(encode)  {}

    /** Number of bytes required  for representing the data
     * \returns size in bytes
     */
    inline size_t getBufferLength() const { return (size() + 3) / 4; }

    /** Affectation operator.
     * \param[in] d : object to be copied.
     * \return the instance
     */
    Data& operator= (const Data& d)
    {
        if (this != &d)
        {
            /** Special case for binary encoding => we have 4 nucleotides in one byte. */
            if (d.getEncoding() == BINARY)
            {
                this->set (d.getBuffer(), d.getBufferLength());
                this->setSize(d.size());
                this->encoding = BINARY;
            }
            else
            {
                this->set (d.getBuffer(), d.size());
                this->encoding = d.getEncoding();
            }
        }
        return *this;
    }

    /** Set the content of this data as a referenced of another Data object.
     * \param[in] ref : referred data
     * \param[in] offset : position to be used in the referred data
     * \param[in] length : length of the data
     */
    void setRef (Data* ref, size_t offset, size_t length)
    {
        /** We call the parent method. */
        Vector<char>::setRef (ref, offset, length);

        /** We set the encoding. */
        encoding = ref->getEncoding();
    }

    /** \copydoc Vector<char>::setRef(char*,size_t) */
    void setRef (char* buffer, size_t length)
    {
        /** We call the parent method. */
        Vector<char>::setRef (buffer, length);
    }

    /** Get the encoding scheme of the data.
     * \return format of the data. */
    Encoding_e getEncoding ()  const  { return encoding; }

    /** Set the encoding scheme of the data.
     * \param[in] encoding : encoding scheme to be used.
     */
    void setEncoding (Encoding_e encoding)  { this->encoding = encoding; }

    /** Conversion from one encoding scheme to another.
     *  TO BE IMPROVED (support only one kind of conversion, from binary to integer)
     * \param[in] in  : input data
     * \param[in] out : output data */
    static void convert (Data& in, Data& out)
    {
        size_t nchar = in.getBufferLength();
        size_t j=0;
        for (size_t i=0; i<nchar; i++)
        {
            char fournt = in[i];
            out[j+3] = fournt & 3; fournt = fournt >> 2;
            out[j+2] = fournt & 3; fournt = fournt >> 2;
            out[j+1] = fournt & 3; fournt = fournt >> 2;
            out[j+0] = fournt & 3;
            j+=4;
        }

        out.encoding = Data::INTEGER;
        out.setSize (in.size());
    }

    /** Shortcut.
     *  - first  : the nucleotide value (A=0, C=1, T=2, G=3)
     *  - second : 0 if valid, 1 if invalid (in case of N character for instance) */
    typedef std::pair<char,char> ConvertChar;

    /** Used to have the following trick:
     *   (buffer[idx]>>3)& 1
     * was equal to 1 for 'A', 'C', 'G' and 'T' (and also the lower case version)
     * and was equal to 0 for 'N' and 'n' 
     * but unfortunately it doesn't work for some of the IUPAC codes, like 'R' 
     * */
    struct ConvertASCII    { static ConvertChar get (const char* buffer, size_t idx)  { return ConvertChar((buffer[idx]>>1) & 3, validNucleotide[(unsigned char)(buffer[idx])]); }};
    struct ConvertInteger  { static ConvertChar get (const char* buffer, size_t idx)  { return ConvertChar(buffer[idx],0); }         };
    struct ConvertBinary   { static ConvertChar get (const char* buffer, size_t idx)  { return ConvertChar(((buffer[idx>>2] >> ((3-(idx&3))*2)) & 3),0); } };

    static const unsigned char validNucleotide[];
private:

    /** Encoding scheme of the data instance. */
    Encoding_e  encoding;

    /* generated using:
     * codes=[1]*256
     * codes[ord('A')] = codes[ord('a')] = 0;
     * codes[ord('C')] = codes[ord('c')] = 0;
     * codes[ord('G')] = codes[ord('g')] = 0;
     * codes[ord('T')] = codes[ord('t')] = 0;
     * print(codes)
     */
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_DATA_HPP_ */
