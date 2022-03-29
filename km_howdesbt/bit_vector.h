#ifndef bit_vector_H
#define bit_vector_H

#include <cstdint>
#include <string>
#include <sdsl/bit_vectors.hpp>
#include <roaring/roaring.h>
#include "bloom_filter_file.h"

// Honor the deployer's choice of custom RRR block size and rank period. This
// choice can be specified on the build command line by defining RRR_BLOCK_SIZE
// and/or RRR_RANK_PERIOD. If no choice has been specified we default to block
// size 255 and rank period 32.
// 
// RRR_BLOCK_SIZE is the number of uncompressed bits in each block of an RRR
// bit vector. This corresponds to the t_bs field of the rrr_vector type
// template in include/sdsl/rrr_vector.hpp. It is typically one of 255, 127,
// 63, etc. We prefer 255, but past versions of sdsl-lite (before April 2017),
// when compiled with clang had a bug that made RRR fail silently with block
// sizes above 127.
//
// To override this default of 255, add -DRRR_BLOCK_SIZE=127 to the Makefile.
//
// RRR_RANK_PERIOD determines how often a rank sample is included in an RRR
// vector. This corresponds to the t_k field of the rrr_vector type template in
// include/sdsl/rrr_vector.hpp. The default is 32, meaning that one rank sample
// is included for every 32nd block. The ability to override this allows us to
// experiment with different rank frequencies and see how performance changes.
// Though the sdsl-lite implementation appears to allow t_k>255, we limit it
// at 255 because we store it in one byte in the bloom filter file.
//
// To override this default of 32, add e.g. -DRRR_RANK_PERIOD=16 to the Makefile.

#ifndef RRR_BLOCK_SIZE
#define RRR_BLOCK_SIZE 255
#endif

#define DEFAULT_RRR_RANK_PERIOD 32
#ifndef RRR_RANK_PERIOD
#define RRR_RANK_PERIOD DEFAULT_RRR_RANK_PERIOD
#endif

#if ((RRR_BLOCK_SIZE<3) || (RRR_BLOCK_SIZE>255))
#error "invalid RRR_BLOCK_SIZE"
#endif

#if ((RRR_RANK_PERIOD<1) || (RRR_RANK_PERIOD>255))
#error "invalid RRR_RANK_PERIOD"
#endif

using sdslbitvector = sdsl::bit_vector;
using sdslrank0     = sdsl::rank_support_v<0>;
using sdslrank1     = sdsl::rank_support_v<1>;
using sdslselect0   = sdsl::select_support_mcl<0>;
using sdslselect1   = sdsl::select_support_mcl<1>;

#define rrr_template_args RRR_BLOCK_SIZE,sdsl::int_vector<>,RRR_RANK_PERIOD
using rrrbitvector  = sdsl::rrr_vector<rrr_template_args>;
using rrrrank0      = sdsl::rank_support_rrr<0,rrr_template_args>;
using rrrrank1      = sdsl::rank_support_rrr<1,rrr_template_args>;
using rrrselect0    = sdsl::select_support_rrr<0,rrr_template_args>;
using rrrselect1    = sdsl::select_support_rrr<1,rrr_template_args>;
#undef rrr_template_args

//using    roarbitvector = roaring_bitmap_t

#define sdslbitvectorHeaderBytes 8	// to grab "raw" bits, skip this many bytes
									// at the start of an sdslbitvector file

//----------
//
// classes in this module--
//
//----------

class BitVector
	{
public:
	BitVector(const std::string& filename, const size_t offset=0, size_t numBytes=0);
	BitVector(const BitVector* srcBv);
	BitVector(std::uint64_t numBits=0);
	virtual ~BitVector();

	virtual std::string class_identity() const { return "BitVector"; }
	virtual std::string identity() const;
	virtual bool modifiable() { return (bits != nullptr); }
	virtual std::uint32_t compressor() const { return bvcomp_uncompressed; }
	virtual bool is_compressed() const { return false; }
	virtual void load();
	virtual void serialized_in(std::ifstream& in);
	virtual void unfinished() {};  // solely for RrrBitVector and RoarBitVector to override
	virtual void finished() {};    // solely for RrrBitVector and RoarBitVector to override
	virtual void save();
	virtual size_t serialized_out(std::ofstream& out, const std::string& filename, const size_t offset=0);
	virtual size_t serialized_out(std::ofstream& out);

	virtual void discard_bits();
	virtual void new_bits(std::uint64_t numBits);
	virtual void replace_bits(sdslbitvector* srcBits);
	virtual void copy_from(const sdslbitvector* srcBits);

	virtual bool is_all_zeros();
	virtual bool is_all_ones();
	virtual void fill(int bitVal);
	virtual void complement();
	virtual void union_with(const sdslbitvector* srcBits);
	virtual void union_with_complement(const sdslbitvector* srcBits);
	virtual void intersect_with(const sdslbitvector* srcBits);
	virtual void mask_with(const sdslbitvector* srcBits);
	virtual void xor_with(const sdslbitvector* srcBits);
	virtual void squeeze_by(const sdslbitvector* srcBits);

	virtual int operator[](std::uint64_t pos) const;
	virtual void write_bit(std::uint64_t pos, int val=1);

	virtual std::uint64_t rank1(std::uint64_t pos);
	virtual std::uint64_t select0(std::uint64_t rank);
	virtual void discard_rank_select();

	virtual std::uint64_t num_bits() const { return numBits; }
	virtual std::uint64_t size() const;

	virtual std::string to_string() const;
	virtual std::string to_complement_string() const;

public:
	bool isResident;
	std::string filename;
	size_t offset;          // offset from start of file to vector's data
	size_t numBytes;        // number of bytes of data; zero means unknown
	sdslbitvector* bits;
	std::uint64_t numBits;  // number of (conceptual) bits in the bit vector
	sdslrank1*   ranker1;
	sdslselect0* selector0;
	std::uint64_t filterInfo; // filter-dependent info for this bit vector;
							// .. typically zero


public:
	static bool       valid_filename (const std::string& filename);
	static std::string compressor_to_string(std::uint32_t compressor);
	static BitVector* bit_vector     (const std::string& filename,
	                                  const std::string& kind="",
	                                  const size_t offset=0,
	                                  const size_t numBytes=0);
	static BitVector* bit_vector     (const std::string& filename,
	                                  const std::uint32_t compressor,
	                                  const size_t offset=0,
	                                  const size_t numBytes=0);
	static BitVector* bit_vector     (const std::uint32_t compressor,
	                                  const std::uint64_t numBits);
	static BitVector* bit_vector     (const std::uint32_t compressor,
	                                  const BitVector* srcBv);
	};


class RrrBitVector: public BitVector
	{
public:
	RrrBitVector(const std::string& filename, const size_t offset=0, size_t numBytes=0, bool readAsUncompressed=false);
	RrrBitVector(const BitVector* srcBv);
	RrrBitVector(std::uint64_t numBits=0);
	virtual ~RrrBitVector();

	virtual std::string class_identity() const { return "RrrBitVector"; }
	virtual std::uint32_t compressor() const { return (writeAsUncompressed)? bvcomp_unc_rrr : bvcomp_rrr; }
	virtual bool is_compressed() const { return (rrrBits != nullptr); }
	virtual void serialized_in(std::ifstream& in);
	virtual void save();
	virtual void unfinished();
	virtual void finished();
	virtual size_t serialized_out(std::ofstream& out);

	virtual void discard_bits();
	virtual void new_bits(std::uint64_t numBits);
	virtual void copy_from(const sdslbitvector* srcBits);
	virtual void copy_from(const rrrbitvector* srcRrrBits);
	virtual void compress();

	virtual bool is_all_zeros();
	virtual bool is_all_ones();

	virtual int operator[](std::uint64_t pos) const;
	virtual void write_bit(std::uint64_t pos, int val=1);

	virtual std::uint64_t rank1(std::uint64_t pos);
	virtual std::uint64_t select0(std::uint64_t rank);
	virtual void discard_rank_select();

	virtual std::uint64_t size() const;

	bool readAsUncompressed;	// true => input file contains the vector still
								// .. in uncompressed form
	bool writeAsUncompressed;	// true => output file is to contain the vector
								// .. still in uncompressed form
	rrrbitvector* rrrBits;		// exclusive of BitVector.bits; at most one of
								// .. bits and rrrBits is non-null at a given
								// .. time
	rrrrank1*   rrrRanker1;		// exclusive of BitVector.ranker1
	rrrselect0* rrrSelector0;	// exclusive of BitVector.selector0
	};


class RoarBitVector: public BitVector
	{
public:
	RoarBitVector(const std::string& filename, const size_t offset=0, size_t numBytes=0, bool readAsUncompressed=false);
	RoarBitVector(const BitVector* srcBv);
	RoarBitVector(std::uint64_t numBits=0);
	virtual ~RoarBitVector();

	virtual std::string class_identity() const { return "RoarBitVector"; }
	virtual std::uint32_t compressor() const { return (writeAsUncompressed)? bvcomp_unc_roar : bvcomp_roar; }
	virtual bool is_compressed() const { return (roarBits != nullptr); }
	virtual void serialized_in(std::ifstream& in);
	virtual void save();
	virtual void unfinished();
	virtual void finished();
	virtual size_t serialized_out(std::ofstream& out);

	virtual void discard_bits();
	virtual void new_bits(std::uint64_t numBits);
	virtual void copy_from(const sdslbitvector* srcBits);
	virtual void copy_from(const roaring_bitmap_t* srcRoarBits);
	virtual void compress();

	virtual bool is_all_zeros();
	virtual bool is_all_ones();

	virtual int operator[](std::uint64_t pos) const;
	virtual void write_bit(std::uint64_t pos, int val=1);

	virtual std::uint64_t rank1(std::uint64_t pos);
	virtual std::uint64_t select0(std::uint64_t rank);
	virtual void discard_rank_select();

	virtual std::uint64_t size() const;

	bool readAsUncompressed;	// true => input file contains the vector still
								// .. in uncompressed form
	bool writeAsUncompressed;	// true => output file is to contain the vector
								// .. still in uncompressed form
	roaring_bitmap_t* roarBits;	// exclusive of BitVector.bits; at most one of
								// .. bits and roarBits is non-null at a given
								// .. time
	};


class RawBitVector: public BitVector
	{
public:
	RawBitVector(const std::string& filename, const size_t offset=0, const std::uint64_t numBits=0);
	// RawBitVector(BitVector* srcBv);  we intentionally omit this constructor
	RawBitVector(std::uint64_t numBits=0);
	virtual ~RawBitVector();

	virtual std::string class_identity() const { return "RawBitVector"; }
	virtual void serialized_in(std::ifstream& in);
	};


class ZerosBitVector: public BitVector
	{
public:
	ZerosBitVector(const std::string& filename, const size_t offset=0, size_t numBytes=0);
	ZerosBitVector(std::uint64_t numBits=0);
	virtual ~ZerosBitVector();

	virtual std::string class_identity() const { return "ZerosBitVector"; }
	virtual std::uint32_t compressor() const { return bvcomp_zeros; }
	virtual bool is_compressed() const { return true; }
	virtual void serialized_in(std::ifstream& in);
	virtual void save();
	virtual size_t serialized_out(std::ofstream& out, const std::string& filename, const size_t offset=0);
	virtual size_t serialized_out(std::ofstream& out);

	virtual void discard_bits();
	virtual void new_bits(std::uint64_t numBits);
	virtual void copy_from(const sdslbitvector* srcBits);
	virtual void copy_from(const rrrbitvector* srcRrrBits);
	virtual void copy_from(const roaring_bitmap_t* srcRoarBits);

	virtual bool is_all_zeros() { return true; }
	virtual bool is_all_ones() { return false; }
	virtual void fill(int bitVal);
	virtual void complement();
	virtual void union_with(const sdslbitvector* srcBits);
	virtual void intersect_with(const sdslbitvector* srcBits);
	virtual void mask_with(const sdslbitvector* srcBits);

	virtual int operator[](std::uint64_t pos) const { return 0; }
	virtual void write_bit(std::uint64_t pos, int val=1);

	virtual std::uint64_t rank1(std::uint64_t pos);
	virtual std::uint64_t select0(std::uint64_t rank);
	virtual void discard_rank_select();

	virtual std::uint64_t size() const;
	};


class OnesBitVector: public ZerosBitVector
	{
public:
	OnesBitVector(const std::string& filename, const size_t offset=0, size_t numBytes=0);
	OnesBitVector(std::uint64_t numBits=0);
	virtual ~OnesBitVector();

	virtual std::string class_identity() const { return "OnesBitVector"; }
	virtual std::uint32_t compressor() const { return bvcomp_ones; }
	virtual bool is_compressed() const { return true; }

	virtual bool is_all_zeros() { return false; }
	virtual bool is_all_ones() { return true; }

	virtual int operator[](std::uint64_t pos) const { return 1; }

	virtual std::uint64_t rank1(std::uint64_t pos);
	virtual std::uint64_t select0(std::uint64_t rank);
	virtual void discard_rank_select();
	};
#endif // bit_vector_H


//----------
//
// support utilities--
//
//----------

void decompress_rrr (const rrrbitvector* rrrBits,
                     void* dstBits, const std::uint64_t numBits);
