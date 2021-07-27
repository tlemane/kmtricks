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

#include <gatb/tools/misc/impl/XmlReader.hpp>

#include <stdarg.h>
#include <stdio.h>

#include <algorithm>
#include <iostream>
#include <string>

#define DEBUG(a)  printf a

using namespace std;

/********************************************************************************/
namespace gatb {  namespace core { namespace tools {  namespace misc {  namespace impl {
/********************************************************************************/

enum State_e { TEXT, OPENING_TAG, TAG_NAME, ATTRIBUTE_NAME, ATTRIBUTE_VALUE };

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
XmlReader::XmlReader (istream& is) : _is(is)
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
XmlReader::~XmlReader ()
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void XmlReader::read ()
{
    processMainState (_is);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void XmlReader::processMainState (std::istream& is)
{
    State_e s  = TEXT;

    string text;

    for (char c=is.get(); is.good(); c=is.get())
    {
        switch (s)
        {
            /************************************************************/
            case TEXT:
            {
                switch (c)
                {
                    case '<':
                    {
                        s = OPENING_TAG;

                        if (text.empty() == false)
                        {
                            /** We first normalize the text. */
                            normalizeText (text);

                            /** We send a notification. */
                            notify (new XmlTagTextEvent (text));

                            /** We clear the text. */
                            text.clear ();
                        }
                        break;
                    }

                    default:
                    {
                        text.append (1, c);
                        break;
                    }
                }
                break;
            }

            /************************************************************/
            case OPENING_TAG:
            {
                is.putback (c);
                processTagState (is);
                s = TEXT;
                break;
            }

            /************************************************************/
            default:
            {
                break;
            }
        };
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void XmlReader::processTagState (std::istream& is)
{
    string tagName;
    string attrName;
    string attrValue;

    bool isClosingTag = false;

    State_e s = TAG_NAME;

    for (char c=is.get(); is.good(); c=is.get())
    {
             if (c == '/')  {  isClosingTag = true;  }
        else if (c == '>')
        {
            if (attrValue.empty() == false)
            {
                notify (new XmlTagAttributeEvent (attrName, attrValue));
                attrValue.clear ();
            }

            if (tagName.empty() == false)
            {
                if (isClosingTag)   { notify (new XmlTagCloseEvent (tagName)); }
                else                { notify (new XmlTagOpenEvent  (tagName)); }
                tagName.clear ();
            }
            break;
        }

        else switch (s)
        {
            /************************************************************/
            case TAG_NAME:
            {
                switch (c)
                {
                    case ' ':
                    {
                        s = ATTRIBUTE_NAME;
                        if (tagName.empty() == false)
                        {
                            notify (new XmlTagOpenEvent (tagName));
                            tagName.clear ();
                        }
                        break;
                    }

                    default:
                    {
                        tagName.append (1, c);
                        break;
                    }
                }
                break;
            }

            /************************************************************/
            case ATTRIBUTE_NAME:
            {
                switch (c)
                {
                    case '=':
                    {
                        s = ATTRIBUTE_VALUE;
                        break;
                    }

                    default:
                    {
                        attrName.append (1, c);
                        break;
                    }
                }
                break;
            }

            /************************************************************/
            case ATTRIBUTE_VALUE:
            {
                switch (c)
                {
                    case ' ':
                    {
                        s = ATTRIBUTE_NAME;
                        if (attrValue.empty() == false)
                        {
                            notify (new XmlTagAttributeEvent (attrName, attrValue));
                            attrName.clear  ();
                            attrValue.clear ();
                        }
                        break;
                    }

                    default:
                    {
                        if (c != '"')  {   attrValue.append (1, c);  }
                        break;
                    }
                }
                break;
            }

            default:
            {
                break;
            }
        }
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void XmlReader::normalizeText (std::string& s)
{
    std::replace (s.begin(), s.end(), '\n', ' ');

    replace (s, "&amp;",   "&");
    replace (s, "&lt;",    "<");
    replace (s, "&gt;",    ">");
    replace (s, "&quot;",  "\"");
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void XmlReader::replace (std::string& str, const std::string& oldStr, const std::string& newStr)
{
  size_t pos = 0;
  while((pos = str.find(oldStr, pos)) != std::string::npos)
  {
     str.replace(pos, oldStr.length(), newStr);
     pos += newStr.length();
  }
}


/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/
