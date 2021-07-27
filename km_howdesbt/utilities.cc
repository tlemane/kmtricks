// utilities.cc-- miscellaneous utility functions.

#include <string>
#include <cstdlib>
#include <cstdint>
#include <climits>
#include <iostream>
#include <vector>
#include <set>

#include "utilities.h"

using std::string;
using std::set;
using std::cerr;
using std::endl;
#define u8  std::uint8_t
#define u32 std::uint32_t
#define u64 std::uint64_t

//----------
//
// nucleotide lookup tables
//
//----------

// ntToComplement maps an *unsigned* ascii character to its complement.
// Upper/lower case are preserved. IUPAC characters (N,SWRYMKBDHV) are also
// supported. Any other characters are unchanged.

static const u8 ntToComplement[256] =
	{
	0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F,
	0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1A, 0x1B, 0x1C, 0x1D, 0x1E, 0x1F,
	0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27, 0x28, 0x29, 0x2A, 0x2B, 0x2C, 0x2D, 0x2E, 0x2F,
	0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3A, 0x3B, 0x3C, 0x3D, 0x3E, 0x3F,
	//     A     B     C     D     E     F     G     H     I     J     K     L     M     N     O
	0x40, 'T',  'V',  'G',  'H',  0x45, 0x46, 'C',  'D',  0x49, 0x4A, 'M',  0x4C, 'K',  'N',  0x4F,
	// P   Q     R     S     T     U     V     W     X     Y     Z
	0x50, 0x51, 'Y',  'S',  'A',  0x55, 'B',  'W',  0x58, 'R',  0x5A, 0x5B, 0x5C, 0x5D, 0x5E, 0x5F,
	//     a     b     c     d     e     f     g     h     i     j     k     l     m     n     o
	0x60, 't',  'v',  'g',  'h',  0x65, 0x66, 'c',  'd',  0x69, 0x6a, 'm',  0x6c, 'k',  'n',  0x6f,
	// p   q     r     s     t     u     v     w     x     y     z
	0x70, 0x71, 'y',  's',  'a',  0x75, 'b',  'w',  0x78, 'r',  0x7a, 0x7b, 0x7c, 0x7d, 0x7e, 0x7f,
	0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89, 0x8A, 0x8B, 0x8C, 0x8D, 0x8E, 0x8F,
	0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98, 0x99, 0x9A, 0x9B, 0x9C, 0x9D, 0x9E, 0x9F,
	0xA0, 0xA1, 0xA2, 0xA3, 0xA4, 0xA5, 0xA6, 0xA7, 0xA8, 0xA9, 0xAA, 0xAB, 0xAC, 0xAD, 0xAE, 0xAF,
	0xB0, 0xB1, 0xB2, 0xB3, 0xB4, 0xB5, 0xB6, 0xB7, 0xB8, 0xB9, 0xBA, 0xBB, 0xBC, 0xBD, 0xBE, 0xBF,
	0xC0, 0xC1, 0xC2, 0xC3, 0xC4, 0xC5, 0xC6, 0xC7, 0xC8, 0xC9, 0xCA, 0xCB, 0xCC, 0xCD, 0xCE, 0xCF,
	0xD0, 0xD1, 0xD2, 0xD3, 0xD4, 0xD5, 0xD6, 0xD7, 0xD8, 0xD9, 0xDA, 0xDB, 0xDC, 0xDD, 0xDE, 0xDF,
	0xE0, 0xE1, 0xE2, 0xE3, 0xE4, 0xE5, 0xE6, 0xE7, 0xE8, 0xE9, 0xEA, 0xEB, 0xEC, 0xED, 0xEE, 0xEF,
	0xF0, 0xF1, 0xF2, 0xF3, 0xF4, 0xF5, 0xF6, 0xF7, 0xF8, 0xF9, 0xFA, 0xFB, 0xFC, 0xFD, 0xFE, 0xFF
	};

// ntIsACGT maps an *unsigned* ascii character to true iff the character is a
// valid A, C, G, or T (upper or lower case).

#define __ false
#define T_ true

static const bool ntIsACGT[256] =
	{
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,
	// A  B  C  D  E  F  G  H  I  J  K  L  M  N  O
	__,T_,__,T_,__,__,__,T_,__,__,__,__,__,__,__,__,
	// Q  R  S  T  U  V  W  X  Y  Z
	__,__,__,__,T_,__,__,__,__,__,__,__,__,__,__,__,
	// a  b  c  d  e  f  g  h  i  j  k  l  m  n  o
	__,T_,__,T_,__,__,__,T_,__,__,__,__,__,__,__,__,
	// q  r  s  t  u  v  w  x  y  z
	__,__,__,__,T_,__,__,__,__,__,__,__,__,__,__,__,
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__
	};

#undef __
#undef T_

//----------
//
// is_prefix_of, is_suffix_of--
//	Determine if one string is a prefix (or suffix) of another.
//
// NOTE: in C++20, callers should use basic_string.starts_with and
//       basic_string.ends_with instead of these functions.
//
//----------
//
// Arguments:
//	const string&	s:		The usually-longer string.
//	const string&	prefix:	the prefix to test for.
//
// Returns:
//	True if s begins with the prefix (or ends with the suffix);  false
//	otherwise. Note that if s is shorter than the prefix it cannot begin with
//	the prefix.
//
//----------

bool is_prefix_of
   (const string&	s,
	const string&	prefix)
	{
	string::size_type prefixLen = prefix.length();
	if (s.length() < prefixLen) return false;
	return (prefix == s.substr(0,prefixLen));
	}

bool is_suffix_of
   (const string&	s,
	const string&	suffix)
	{
	string::size_type sLen      = s.length();
	string::size_type suffixLen = suffix.length();
	if (sLen < suffixLen) return false;
	return (suffix == s.substr(sLen-suffixLen));
	}

//----------
//
// strip_blank_ends, strip_blank_prefix, strip_blank_suffix--
//	Remove the blank prefix and/or suffix from a string.
//
//----------
//
// Arguments:
//	const string&	s:	The string to examine.
//
// Returns:
//	A copy of s with the blank end(s) removed.
//
//----------

string strip_blank_ends
   (const string& s)
	{
	return strip_blank_prefix(strip_blank_suffix(s));
	}

string strip_blank_prefix
   (const string& s)
	{
	string ss = s;
	return ss.erase(0, ss.find_first_not_of(" "));
	}

string strip_blank_suffix
   (const string& s)
	{
	string ss = s;
	return ss.erase(ss.find_last_not_of(" ") + 1);
	}

//----------
//
// strip_prefix, strip_suffix--
//	Remove a specified prefix or suffix from a string.
//
//----------
//
// Arguments:
//	const string&	s:	The string to examine.
//
// Returns:
//	A copy of s with the specified prefix or suffix removed. If the string did
//	not contain the prefix or suffix, the original string is returned.
//
//----------

string strip_prefix
   (const string& s,
	const string& prefix)
	{
	if (is_prefix_of (s, prefix))
		return s.substr(prefix.length());
	else
		return s;
	}

string strip_suffix
   (const string& s,
	const string& suffix)
	{
	if (is_suffix_of (s, suffix))
		return s.substr(0,s.length()-suffix.length());
	else
		return s;
	}

//----------
//
// strip_file_path--
//	Remove any path prefix from a filename.
//
//----------
//
// Arguments:
//	const string&	filename:	The string to examine.
//
// Returns:
//	A copy of filename with any path prefix removed.
//
//----------

string strip_file_path
   (const string& filename)
	{
	string::size_type slashIx = filename.find_last_of("/");
	if (slashIx == string::npos)
		return filename;
	else
		return filename.substr(slashIx+1);
	}

//----------
//
// string_to_int, string_to_u32, string_to_u64--
//	Parse a string for the integer value it contains.
//
//----------
//
// Arguments:
//	const string&	s:			The string to parse.
//	const bool		allowHex:	(true) parse strings beginning with "0x" as
//								.. hexadecimal.
//
// Returns:
//	The integer value of the string. Note that the string *must not* contain
//	anything other than a valid integer-- failures result in program
//	termination.
//
//----------

int string_to_int
   (const string&	s,
	const bool		allowHex)
	{
	int		v;
	size_t	extra;

	// an empty string is not a number
	// $$$ we should check whether the string just contains blanks

	if (s == "") goto empty_string;

	// if it begins with 0x, parse it as hex
	// nota bene: we avoid INT_MIN because on some platforms, it is
	//            apparently defined as unsigned(!)

	if ((allowHex) && (is_prefix_of (s, "0x")))
		{
		std::int64_t vv = (std::int64_t) hex_string_to_u64 (s);
		if ((vv >= 0) && ((u64)  vv  <=  (u64) INT_MAX))    return (int) vv;
		if ((vv <  0) && ((u64)(-vv) <= ((u64) INT_MAX+1))) return (int) vv;
		goto out_of_range;
		}

	// use stoi to parse it, but make sure stoi used the entire string
	// $$$ we should check whether the leftover suffix is just blanks

	v = std::stoi (s, &extra);
	if (extra != s.length()) goto not_an_integer;

	return v;

	//////////
	// failure exits
	//////////

empty_string:
	fatal ("an empty string is not an integer");

not_an_integer:
	fatal ("\"" + s + "\" is not an integer");

out_of_range:
	fatal ("\"" + s + "\" is out of range for an integer");
	return 0;  // execution never reaches here
	}


u32 string_to_u32
   (const string&	s,
	const bool		allowHex)
	{
	u64 v = string_to_u64 (s, allowHex);

	if (v <= (u64) UINT32_MAX)
		return (u32) v;

	fatal ("\"" + s + "\" is out of range for a 32-bit unsigned integer");
	return 0;  // execution never reaches here
	}


u64 string_to_u64
   (const string&	s,
	const bool		allowHex)
	{
	unsigned long long	v;
	size_t				extra;

	// an empty string is not a number
	// $$$ we should check whether the string just contains blanks

	if (s == "") goto empty_string;

	// if it begins with 0x, parse it as hex

	if ((allowHex) && (is_prefix_of (s, "0x")))
		return hex_string_to_u64 (s);

	// use stoull to parse it, but make sure stoull used the entire string
	// $$$ we should check whether the leftover suffix is just blanks

	v = std::stoull (s, &extra);
	if (extra != s.length()) goto not_an_integer;

	return v;

	//////////
	// failure exits
	//////////

empty_string:
	fatal ("an empty string is not an unsigned integer");

not_an_integer:
	fatal ("\"" + s + "\" is not an unsigned integer");
	return 0;  // execution never reaches here
	}

//----------
//
// string_to_unitized_int, string_to_unitized_u32, string_to_unitized_u64--
//	Parse a string for the integer value it contains, allowing a unit suffix,
//	e.g. K, M, G, or T.
//
//----------
//
// Arguments:
//	const string&	s:			The string to parse.
//	const int		unitScale:	The unit multiplier, which must be either 1000
//								or 1024.
//
// Returns:
//	The integer value of the string. Note that the string *must not* contain
//	anything other than a valid integer (with or without a unit suffix)-- failures result in program
//	termination.
//
//----------

int string_to_unitized_int
   (const string&	s,
	const int		_unitScale)
	{
	string			parseMe = s;
	size_t			sLen = s.length();
	u64				unitScale = _unitScale;
	u64				multiplier;
	bool			isFloat;
	int				v;
	double			vf = 0.0;
	size_t			extra;

	if ((_unitScale != 1000) && (_unitScale != 1024)) goto bad_unit_scale;

	// an empty string is not a number
	// $$$ we should check whether the string just contains blanks

	if (s == "") goto empty_string;

	// extract any unit from the end

	multiplier = 1;
	switch (s[sLen-1])
		{
		case 'G': case 'g': multiplier *= unitScale;
		case 'M': case 'm': multiplier *= unitScale;
		case 'K': case 'k': multiplier *= unitScale;
			break;
		default:
			break;
		}

	if (multiplier != 1)
		{
		if (sLen == 1) goto not_an_integer; // (empty except for the unit)
		parseMe = s.substr (0, sLen-1);
		sLen--;
		}

	// use stoi to parse it first, but if that fails to consume the entire
	// string, try parsing it as a float with stod (but only if we had a unit)

	isFloat = false;
	v = std::stoi (parseMe, &extra);
	if ((multiplier != 1) && (extra != sLen))
		{
		vf = std::stod (parseMe, &extra);
		if (extra != sLen) goto not_an_integer;
		isFloat = true;
		}

	// apply the multiplier, being careful not to overflow

	if (isFloat)
		{
		if ((vf > 0) && (vf*multiplier > INT32_MAX)) goto overflow;
		if ((vf < 0) && (vf*multiplier < INT32_MIN)) goto overflow;
		v = (vf * multiplier) + .5;
		}
	else if (multiplier != 1)
		{
		// nota bene: we avoid INT32_MIN because on some platforms, it is
		//            apparently defined as unsigned(!)
		if ((v > 0) && (((u32)  v) >  ((u32)INT32_MAX)    / multiplier)) goto overflow;
		if ((v < 0) && (((u32) -v) > (((u32)INT32_MAX)+1) / multiplier)) goto overflow;
		v *= multiplier;
		}

	return v;

	//////////
	// failure exits
	//////////

bad_unit_scale:
	fatal ("internal error: string_to_unitized_int(*," + std::to_string(_unitScale) + ")");

empty_string:
	fatal ("an empty string is not an integer");

not_an_integer:
	fatal ("\"" + s + "\" is not an integer");

overflow:
	fatal ("\"" + s + "\" is out of range for an integer");
	return 0;  // execution never reaches here
	}


u32 string_to_unitized_u32
   (const string&		s,
	const int			unitScale)
	{
	if ((unitScale != 1000) && (unitScale != 1024))
		fatal ("internal error: string_to_unitized_u32(*," + std::to_string(unitScale) + ")");

	u64 v = string_to_unitized_u64 (s, unitScale);

	if (v <= (u64) UINT32_MAX)
		return (u32) v;

	fatal ("\"" + s + "\" is out of range for a 32-bit unsigned integer");
	return 0;  // execution never reaches here
	}


u64 string_to_unitized_u64
   (const string&		s,
	const int			_unitScale)
	{
	string				parseMe = s;
	size_t				sLen = s.length();
	u64					unitScale = _unitScale;
	u64					multiplier;
	bool				isFloat;
	unsigned long long	v;
	double				vf = 0.0;
	size_t				extra;

	if ((_unitScale != 1000) && (_unitScale != 1024)) goto bad_unit_scale;

	// an empty string is not a number
	// $$$ we should check whether the string just contains blanks

	if (s == "") goto empty_string;

	// extract any unit from the end
	// nota bene: we don't implement E for exa here, because E is often used
	//            as a prefix in scientific notation, e.g. "1e-6"; though this
	//            wouldn't cause any real parsing difficulty, we don't expect
	//            things like exabytes are in common use

	multiplier = 1;
	switch (s[sLen-1])
		{
		case 'P': case 'p': multiplier *= unitScale;
		case 'T': case 't': multiplier *= unitScale;
		case 'G': case 'g': multiplier *= unitScale;
		case 'M': case 'm': multiplier *= unitScale;
		case 'K': case 'k': multiplier *= unitScale;
			break;
		default:
			break;
		}

	if (multiplier != 1)
		{
		if (sLen == 1) goto not_an_integer; // (empty except for the unit)
		parseMe = s.substr (0, sLen-1);
		sLen--;
		}

	// use stoull to parse it first, but if that fails to consume the entire
	// string, try parsing it as a float with stod (but only if we had a unit)

	isFloat = false;
	v = std::stoull (parseMe, &extra);

	if (extra != sLen)
		{
		if (multiplier == 1) goto not_an_integer;
		vf = std::stod (parseMe, &extra);
		if (extra != sLen) goto not_an_integer;
		isFloat = true;
		}

	// apply the multiplier, being careful not to overflow

	if (isFloat)
		{
		if (vf < 0) goto overflow;
		if (vf * multiplier > UINT64_MAX) goto overflow;
		v = (vf * multiplier) + .5;
		}
	else if (multiplier != 1)
		{
		if (v > UINT64_MAX / multiplier) goto overflow;
		v *= multiplier;
		}

	return v;

	//////////
	// failure exits
	//////////

bad_unit_scale:
	fatal ("internal error: string_to_unitized_u64(*," + std::to_string(_unitScale) + ")");

empty_string:
	fatal ("an empty string is not an unsigned integer");

not_an_integer:
	fatal ("\"" + s + "\" is not an unsigned integer");

overflow:
	fatal ("\"" + s + "\" is out of range for an unsigned integer");
	return 0;  // execution never reaches here
	}

//----------
//
// hex_string_to_u32, hex_string_to_u64--
//	Parse a string for the hexadecimal integer value it contains.
//
//----------
//
// Arguments:
//	const string&	s:	The string to parse. This may have an optional "0x"
//						prefix.
//
// Returns:
//	The integer value of the string. Note that the string *must not* contain
//	anything other than a valid hexadecimal integer-- failures result in
//	program termination.
//
//----------

u32 hex_string_to_u32
   (const string&	s)
	{
	unsigned long	v;
	size_t			extra;

	// strip the optional "0x" prefix

	string parseMe = s;
	if (is_prefix_of (parseMe, "0x"))
		parseMe = parseMe.substr(2);

	// an empty string is not a number
	// $$$ we should check whether the string just contains blanks

	if (parseMe == "") goto empty_string;

	// use stoul to parse it, but make sure stoul used the entire string
	// $$$ we should check whether the leftover suffix is just blanks

	v = std::stoul (parseMe, &extra, 16);
	if (extra != parseMe.length()) goto not_an_integer;

	return v;

	//////////
	// failure exits
	//////////

empty_string:
	fatal ("an empty string is not an hexadecimal unsigned integer");

not_an_integer:
	fatal ("\"" + s + "\" is not an hexadecimal unsigned integer");
	return 0;  // execution never reaches here
	}

u64 hex_string_to_u64
   (const string&	s)
	{
	unsigned long	v;
	size_t			extra;

	// strip the optional "0x" prefix

	string parseMe = s;
	if (is_prefix_of (parseMe, "0x"))
		parseMe = parseMe.substr(2);

	// an empty string is not a number
	// $$$ we should check whether the string just contains blanks

	if (parseMe == "") goto empty_string;

	// use stoull to parse it, but make sure stoull used the entire string
	// $$$ we should check whether the leftover suffix is just blanks

	v = std::stoull (parseMe, &extra, 16);
	if (extra != parseMe.length()) goto not_an_integer;

	return v;

	//////////
	// failure exits
	//////////

empty_string:
	fatal ("an empty string is not an hexadecimal unsigned integer");

not_an_integer:
	fatal ("\"" + s + "\" is not an hexadecimal unsigned integer");
	return 0;  // execution never reaches here
	}

//----------
//
// string_to_double--
//	Parse a string for the floating point number it contains.
//
// Values can be expressed as real numbers, percentages, or fractions. Examples
// are "0.3", "30%" and "3/10".
//
//----------
//
// Arguments:
//	const string&	s:	The string to parse.
//
// Returns:
//	The value of the string. Note that the string *must not* contain anything
//	other than a valid floating point number-- failures result in program
//	termination.
//
//----------

double string_to_double
   (const string&	s)
	{
	double v;
	size_t extra;
	string leftover;

	// an empty string is not a number
	// $$$ we should check whether the string just contains blanks

	if (s == "") goto empty_string;

	// use stod to parse it, but make sure stod used the entire string
	// $$$ we should check whether the leftover suffix is just blanks

	v = std::stod (s, &extra);
	if (extra == s.length()) return v;

	leftover = s.substr(extra);

	// if the leftover part is a %, treat the value as a percentage

	if (leftover == "%")
		{ v /= 100.0;  return v; }

	// if the leftover part is a /, treat the value as a numerator, and parse
	// the denominator
	// $$$ we shouldn't allow 0 as a denominator

	if (leftover[0] == '/') // (nota bene: leftover cannot be empty)
		{
		leftover = leftover.substr(1);
		if (leftover == "") goto not_a_number;
		double denom = std::stod (leftover, &extra);
		if (extra != leftover.length()) goto not_a_number;
		v /= denom;
		return v;
		}

	// otherwise, it's not a valid number

	goto not_a_number;

	//////////
	// failure exits
	//////////

empty_string:
	fatal ("an empty string is not a number");

not_a_number:
	fatal ("\"" + s + "\" is not a valid number");
	return 0;  // execution never reaches here
	}

//----------
//
// string_to_probability--
//	Parse a string for the probability value it contains.
//
// Probabilities can be expressed as real numbers, percentages, or fractions,
// but must be in the interval 0..1. Examples are "0.3", "30%" and "3/10".
//
//----------
//
// Arguments:
//	const string&	s:	The string to parse.
//
// Returns:
//	The value of the string. Note that the string *must not* contain anything
//	other than a valid probability-- failures result in program termination.
//
//----------

double string_to_probability
   (const string&	s)
	{
	double v;
	size_t extra;
	string leftover;

	// an empty string is not a number
	// $$$ we should check whether the string just contains blanks

	if (s == "") goto empty_string;

	// use stod to parse it, but make sure stod used the entire string
	// $$$ we should check whether the leftover suffix is just blanks

	v = std::stod (s, &extra);
	if (extra == s.length()) goto verify_probability;

	leftover = s.substr(extra);

	// if the leftover part is a %, treat the value as a percentage

	if (leftover == "%")
		{ v /= 100.0;  goto verify_probability; }

	// if the leftover part is a /, treat the value as a numerator, and parse
	// the denominator
	// $$$ we shouldn't allow 0 as a denominator

	if (leftover[0] == '/') // (nota bene: leftover cannot be empty)
		{
		leftover = leftover.substr(1);
		if (leftover == "") goto not_a_probability;
		double denom = std::stod (leftover, &extra);
		if (extra != leftover.length()) goto not_a_probability;
		v /= denom;
		goto verify_probability;
		}

	// otherwise, it's not a valid probability

	goto not_a_probability;

	// make sure the probability is in the unit interval

verify_probability:
	if ((v < 0.0) || (v > 1.0)) goto not_a_probability;
	return v;

	//////////
	// failure exits
	//////////

empty_string:
	fatal ("an empty string is not a probability");

not_a_probability:
	fatal ("\"" + s + "\" is not a valid probability");
	return 0;  // execution never reaches here
	}

//----------
//
// to_lower--
//	Create a lowercase copy of a string.
//
//----------
//
// Arguments:
//	const string&	s:	The string to copy.
//
// Returns:
//	A copy of the string, in lowercase.
//
//----------

string to_lower
   (const string&	s)
	{
	string::size_type sLen = s.length();
	string	lc(sLen,' ');

	string::size_type ix, iy;
	for (ix=0,iy=sLen-1 ; ix<sLen ; ix++,iy--)
		lc[ix] = tolower(s[ix]);

	return lc;
	}

//----------
//
// reverse_complement--
//	Create the reverse complement of a nucleotide string.
//
//----------
//
// Arguments:
//	const string&	s:	The string to reverse complement.
//
// Returns:
//	The reverse complement of the string.
//
//----------

string reverse_complement
   (const string&	s)
	{
	string::size_type sLen = s.length();
	string	rc(sLen,' ');

	string::size_type ix, iy;
	for (ix=0,iy=sLen-1 ; ix<sLen ; ix++,iy--)
		rc[ix] = (char) ntToComplement[(u8)s[iy]];

	return rc;
	}

//----------
//
// nt_is_acgt--
//	Check whether an alleged nucleotide is really a nucleotide.
//
//----------
//
// Arguments:
//	const char nt:	The (alleged) nucleotide to check.
//
// Returns:
//	true if the nt is one of A, C, G, T (upper or lower case); false otherwise.
//
//----------

bool nt_is_acgt
   (const char nt)
	{
	return ntIsACGT[(u8)nt];
	}

//----------
//
// update_crc--
//	Incorporate the next byte into a cyclic redundancy check.
//
//----------
//
// Arguments:
//	u32 crc:	The crc of the previous bytes
//	u8	ch:		The byte to "add" to the crc.
//
// Returns:
//	The new crc value.
//
//----------
//
// Notes:
//	(1)	This code is based on some borrowed ages from
//		  http://remus.rutgers.edu/~rhoads/Code/crc-32b.c
//	(2)	To compute a crc over a string of bytes, do something like this:
//		  crc = 0
//		  for ch in string:
//		    crc = update_crc(crc,ch)
//
//----------

static bool crcTableInitialized = false;
static u32 crcTable[256];


static void generate_crc_table (void)
	{
	u32	crc, poly;
	int	i, j;

	if (crcTableInitialized) return;
	crcTableInitialized = true;
	
	poly = 0xEDB88320L;
	for (i=0 ; i<256 ; i++)
		{
		crc = i;
		for (j=8 ; j>0 ; j--)
			{
			if (crc & 1) crc = (crc >> 1) ^ poly;
			        else crc >>= 1;
			}
		crcTable[i] = crc;
		}
	}


u32 update_crc (u32 crc, u8 ch)
	{
	if (not crcTableInitialized) generate_crc_table();
	return (crc>>8) ^ crcTable[(crc^ch) & 0xFF];
	}

//----------
//
// fatal--
//	Cause program fatality, after pushing a message out to the user.
//
//----------

void fatal
   (const string& message)
	{
	if (!message.empty()) cerr << message << endl;
	std::exit (EXIT_FAILURE);
	}
