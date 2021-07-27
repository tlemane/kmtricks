// support.cc-- miscellaneous support functions.

#include <string>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include "support.h"

using std::string;
using std::vector;

//----------
//
// parse_comma_list--
//	Parse a string for the list of comma-separated fields it contains.
//
//----------
//
// Arguments:
//	const string&	s:	The string to parse.
//
// Returns:
//	A vector of the fields in the string.
//
//----------

vector<string> parse_comma_list
   (const string&	s)
	{
	vector<string>	fields;

	std::istringstream ss(s);
	while (ss)
		{
		string field;
		if (!getline (ss, field, ',' )) break;
		fields.emplace_back (field);
		}

	return fields;
	}

//----------
//
// tokenize, quoted_tokenize--
//	Break a string into its whitespace-separated fields.
//
//----------
//
// Arguments:
//	const string&	s:		The string to reverse complement.
//
// Returns:
//	A vector of the tokens in the string.
//
//----------
//
// Notes:
//	(1)	quoted_tokenize() recognizes quoted strings, and is more suitable for
//		command-line parsing.
//	(2)	tokenize() is based on one of the solutions posted at
//		  http://stackoverflow.com/questions/53849/how-do-i-tokenize-a-string-in-c
//
//----------

vector<string> tokenize
   (const string&	s)
	{
	vector<string>	tokens;
	string::size_type startIx, endIx;

	startIx = s.find_first_not_of (" \t\n");
	endIx   = startIx;

	while (startIx != std::string::npos)
		{
		endIx = s.find_first_of (" \t\n", startIx);
		if (endIx == std::string::npos) endIx = s.length();
		tokens.emplace_back (s.substr (startIx, endIx-startIx));
		startIx = s.find_first_not_of (" \t\n", endIx);
		}

	return tokens;
	}


vector<string> quoted_tokenize
   (const string&	s)
	{
	enum { whitespace = 0, darkspace, quoted } state;
	vector<string>	tokens;

	char token[s.length()+1];
	string::size_type tokenLen = 0;

	state = whitespace;
	for (string::size_type ix=0 ; ix<s.length() ; ix++)
		{
		char ch = s[ix];
		bool finishToken = false;
		switch (state)
			{
			default:
			case whitespace:
				if (ch == ' ')    break;
				if (ch == '\t')   break;
				if (ch == '\"') { state = quoted;                              break; }
				if (ch != '\\') { state = darkspace;  token[tokenLen++] = ch;  break; }

				if (ix+1 >= s.length()) break;  // escape at end of string is ignored
				ch = s[ix++];
				state = darkspace;  token[tokenLen++] = ch;
				break;

			case darkspace:
				if (ch == ' ')  { finishToken = true;      break; }
				if (ch == '\t') { finishToken = true;      break; }
				if (ch == '\"') { state = quoted;          break; }
				if (ch != '\\') { token[tokenLen++] = ch;  break; }

				if (ix+1 >= s.length()) break;  // escape at end of string is ignored
				ch = s[ix++];
				token[tokenLen++] = ch;
				break;

			case quoted:
				if (ch == '\"') { state = darkspace;       break; }
				if (ch != '\\') { token[tokenLen++] = ch;  break; }

				if (ix+1 >= s.length()) break;  // escape at end of string is ignored
				ch = s[ix++];
				token[tokenLen++] = ch;
				break;
			}

		if (finishToken)
			{
			if (tokenLen > 0)
				{
				token[tokenLen] = 0;
				tokens.emplace_back (string(token));
				}
			state = whitespace;  // nota bene: andy time we switch to whitespace,
			tokenLen = 0;        // .. we have to set tokenLen to zero
			}
		}

	if (tokenLen > 0)
		{
		token[tokenLen] = 0;
		tokens.emplace_back (string(token));
		}

	return tokens;
	}

//----------
//
// expand_filenames--
//	Copy a list of filenames. Names containing {number} are copied multiple
//	times, replacing {number} with a different number for each copy (from 1 to
//	fileCount).
//
//----------
//
// Arguments:
//	const vector<string>	filenames:	The filenames to copy.
//	int						fileCount:	How many copies to make of any
//										.. filenames containing {number}.
//	vector<string>&			filenames:	List to copy filenames to.
//
// Returns:
//	(nothing)
//
//----------

void expand_filenames
   (const vector<string>	filenames,
	int						fileCount,
	vector<string>&			newFilenames)
	{
	for (const auto& filename : filenames)
		{
		string field = "{number}";
		std::size_t fieldIx = filename.find (field);
		if (fieldIx == string::npos)
			newFilenames.emplace_back(filename);
		else
			{
			for (int fileNum=1 ; fileNum<=fileCount ; fileNum++)
				{
				string name = filename;
				name.replace (fieldIx, field.length(), std::to_string(fileNum));
				newFilenames.emplace_back(name);
				}
			}
		}

	}
