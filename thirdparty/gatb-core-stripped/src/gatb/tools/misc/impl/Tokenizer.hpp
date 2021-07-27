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

/** \file Tokenizer.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief String Tokenizer
 */

#ifndef _GATB_CORE_TOOLS_MISC_TOKENIZER_HPP_
#define _GATB_CORE_TOOLS_MISC_TOKENIZER_HPP_

/********************************************************************************/

#include <gatb/system/api/types.hpp>
#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/** \brief  String tokenizer as an Iterator.
 *
 *  Tool for tokenizing strings (like strtok) that follows our Iterator concept.
 *
 *  One should provide both the string to be tokenized and a string holding
 *  characters to be considered as separators.
 *
 *  Code sample:
 *  \code
 *  void sample ()
 *  {
 *      // We create a tokenizer with spaces as separators.
 *      TokenizerIterator it ("this is the text to tokenize", " ");
 *
 *      // We loop the iterator
 *      for (it.first(); !it.isDone(); it.next())
 *      {
 *          // We get the current token
 *          char* token = it.currentItem ();
 *      }
 *  }
 *  \endcode
 */
class TokenizerIterator : public tools::dp::Iterator<char*>
{
public:

    /** Constructors.
     * \param text : the string to be tokenized.
     * \param separator : the separator characters.
     */
    TokenizerIterator (const char* text, const char* separator);

    /** Destructor. */
    virtual ~TokenizerIterator ();

    /** \copydoc tools::dp::Iterator::first */
    void first();

    /** \copydoc tools::dp::Iterator::next */
    void next();

    /** \copydoc tools::dp::Iterator::isDone */
    bool isDone()  { return _loop == 0;  }

    /** \copydoc tools::dp::Iterator::item */
    char*& item()  {  return _loop;  }

private:

    /** String holding separator characters. */
    std::string _sep;

    /** Internal string. */
    char* _text;

    /** Internal string. */
    char* _str;

    /** Internal string. */
    char* _loop;

    /** Internal string. */
    char* _save;

    /** */
    char* tok_r (char* s, const char* delim, char** lasts);
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_TOKENIZER_HPP_ */
