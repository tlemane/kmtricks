#ifndef utilities_H
#define utilities_H

#include <string>
#include <vector>
#include <set>
#include <algorithm>

//----------
//
// prototypes for functions in this module--
//
//----------

bool          is_prefix_of           (const std::string& s, const std::string& prefix);
bool          is_suffix_of           (const std::string& s, const std::string& suffix);
std::string   strip_blank_ends       (const std::string& s);
std::string   strip_blank_prefix     (const std::string& s);
std::string   strip_blank_suffix     (const std::string& s);
std::string   strip_prefix           (const std::string& s, const std::string& prefix);
std::string   strip_suffix           (const std::string& s, const std::string& suffix);
std::string   strip_file_path        (const std::string& filename);
int           string_to_int          (const std::string& s, const bool allowHex=false);
std::uint32_t string_to_u32          (const std::string& s, const bool allowHex=false);
std::uint64_t string_to_u64          (const std::string& s, const bool allowHex=false);
int           string_to_unitized_int (const std::string& s, const int unitScale=1000);
std::uint32_t string_to_unitized_u32 (const std::string& s, const int unitScale=1000);
std::uint64_t string_to_unitized_u64 (const std::string& s, const int unitScale=1000);
std::uint32_t hex_string_to_u32      (const std::string& s);
std::uint64_t hex_string_to_u64      (const std::string& s);
double        string_to_double       (const std::string& s);
double        string_to_probability  (const std::string& s);
std::string   to_lower               (const std::string& s);
std::string   reverse_complement     (const std::string& s);
bool          nt_is_acgt             (const char nt);
std::uint32_t update_crc             (std::uint32_t crc, std::uint8_t ch);
void          fatal                  (const std::string& message="");

#define round_up_16(b)  ((((std::uint64_t) (b))+15)&(~15))

// timer stuff

#include <chrono>
#define get_wall_time()          std::chrono::system_clock::now()
#define elapsed_wall_time(start) ((std::chrono::duration_cast<std::chrono::nanoseconds>(get_wall_time()-start).count()) / 1000000000.0)
#define wall_time_ty             std::chrono::time_point<std::chrono::system_clock>


// macro to declare a variable of type T with the gnu c compiler not caring
// that it's not used

#ifdef __GNUC__
#define notused(T) T unused##__LINE__ __attribute__ ((unused))
#else
#define notused(T) T unused##__LINE__
#endif // __GNUC__

// macro to convince gnu c compiler not to complain about unusued function
// arguments

#ifdef __GNUC__
#define arg_dont_complain(arg) arg __attribute__ ((unused))
#else
#define arg_dont_complain(arg) arg
#endif // __GNUC__

//----------
//
// contains--
//	Determine if a set contains a particular element.
//
//----------
//
// Arguments (variant 1):
//	const std::set<T>&		container:	The set.
//	const T&				element:	The element to check for.
//
// Arguments (variant 2):
//	const std::vector<T>&	container:	The set, as a vector.
//	const T&				element:	The element to check for.
//
// Returns:
//	true if the set (or vector) contains the element;  false otherwise.
//
//----------
//
// Notes:
//	[1]	It's hard to understand why C++ doesn't provide this function.
//	[2]	This code was based on the discussions at
//		  http://stackoverflow.com/questions/1701067/how-to-check-that-an-element-is-in-a-stdset
//		and
//		  http://stackoverflow.com/questions/6277646/in-c-check-if-stdvectorstring-contains-a-certain-value
//
//----------

template<class T>
bool contains
   (const std::set<T>&	container,
	const T&			element)
	{
	return (container.find(element) != container.end());
	}

template<class T>
bool contains
   (const std::set<T>&	container,
	const char*			element)
	{
	return (container.find(std::string(element)) != container.end());
	}

template<class T>
bool contains
   (const std::vector<T>&	container,
	const T&				element)
	{
	return (std::find(container.begin(),container.end(),element) != container.end());
	}

template<class T>
bool contains
   (const std::vector<T>&	container,
	const char*				element)
	{
	return (std::find(container.begin(),container.end(),std::string(element)) != container.end());
	}

#endif // utilities_H
