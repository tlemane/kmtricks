#ifndef bloom_filter_H
#define bloom_filter_H

#include <cstdint>
#include <string>
#include <iostream>
#include <vector>
#include "hash.h"
#include "bit_vector.h"

class FileManager;

//----------
//
// classes in this module--
//
//----------

class BloomFilter
	{
public:
	static const int maxBitVectors = 2;

public:
	BloomFilter(const std::string& filename);
	BloomFilter(const std::string& filename, uint32_t kmerSize,
	            std::uint32_t numHashes, std::uint64_t hashSeed1, std::uint64_t hashSeed2,
	            std::uint64_t numBits, std::uint64_t hashModulus=0);
	BloomFilter(const BloomFilter* templateBf, const std::string& newFilename="");
	virtual ~BloomFilter();

	virtual std::string class_identity() const { return "BloomFilter"; }
	virtual std::string identity() const;
	virtual std::uint32_t kind() const { return bfkind_simple; }

	virtual void setup_hashers();
	virtual bool preload(bool bypassManager=false,bool stopOnMultipleContent=true);
	virtual void load(bool bypassManager=false,const std::string& whichNodeName="");
	virtual void save();

	virtual void copy_properties(const BloomFilter* templateBf);
	virtual void steal_bits(BloomFilter* templateBf);
	virtual void steal_bits(BloomFilter* templateBf, int whichBv);
	virtual void steal_bits(BloomFilter* templateBf, int whichSrcBv, int whichDstBv, std::uint32_t compressor=bvcomp_uncompressed);

	virtual bool is_consistent_with (const BloomFilter* bf, bool beFatal=false) const;

	virtual void discard_bits();
	virtual void discard_bits(int whichBv);
	virtual void new_bits (std::uint32_t compressor=bvcomp_uncompressed, int whichBv=-1);
	virtual void new_bits (BitVector* srcBv, std::uint32_t compressor=bvcomp_uncompressed, int whichBv=0);
	virtual void new_bits (const std::string& filename);
	virtual BitVector* get_bit_vector(int whichBv=0) const;
	virtual BitVector* surrender_bit_vector(int whichBv);
	virtual BitVector* simplify_bit_vector(int whichBv);

	virtual void complement (int whichDstBv=0);
	virtual void union_with (BitVector* srcBv, int whichDstBv=0);
	virtual void union_with_complement (BitVector* srcBv, int whichDstBv=0);
	virtual void intersect_with (BitVector* srcBv, int whichDstBv=0);
	virtual void intersect_with_complement (BitVector* srcBv, int whichDstBv=0) { return mask_with(srcBv,whichDstBv); }
	virtual void mask_with (BitVector* srcBv, int whichDstBv=0);
	virtual void xor_with (BitVector* srcBv, int whichDstBv=0);
	virtual void squeeze_by (BitVector* srcBv, int whichDstBv=0);
	virtual void squeeze_by (const sdslbitvector* srcBits, int whichDstBv=0);

	virtual std::uint64_t mer_to_position(const std::string& mer) const;
	virtual std::uint64_t mer_to_position(const std::uint64_t* merData) const;
	virtual void add (const std::string& mer);
	virtual void add (const std::uint64_t* merData);
	virtual bool contains (const std::string& mer) const;
	virtual bool contains (const std::uint64_t* merData) const;
	virtual int lookup (const std::uint64_t pos) const;

	virtual std::uint64_t hash_modulus() const { return hashModulus; }
	virtual std::uint64_t num_bits()     const { return numBits; }

	virtual bool is_position_adjustor () { return false; }
	virtual void adjust_positions_in_list  (std::vector<std::uint64_t> &kmerPositions,
	                                        std::uint64_t numUnresolved) {}
	virtual void restore_positions_in_list (std::vector<std::uint64_t> &kmerPositions,
	                                        std::uint64_t numUnresolved) {}

public:
	bool ready;						// ready is false until we know the bloom
									// filter's attributes (e.g. kmerSize, hash
									// functions, etc.)
	FileManager* manager;
	std::string filename;
	std::uint32_t kmerSize;

	HashCanonical* hasher1, *hasher2;
	std::uint32_t numHashes;		// how many hash values are generated for
									// .. each key
	std::uint64_t hashSeed1, hashSeed2;

	std::uint64_t hashModulus;		// hash output is reduced to 0..hashModulus-1
	std::uint64_t numBits;			// how many hashed positions are populated
									// .. in each bit vector
									// .. hashModulus >= numBits >= 2

	bool          setSizeKnown;		// true  => the setSize field is valid
									// false => the value of setSize is unknown
	std::uint64_t setSize;			// number of distinct kmers that were
									// .. inserted during construction

	int numBitVectors;				// how many bit vectors are in use
									// .. nota bene: in this default class,
									// .. numBitVectors is 1
	BitVector* bvs[maxBitVectors];

public:
	bool reportLoad = false;
	bool reportSave = false;
	static bool reportSimplify;
	static bool reportLoadTime;
	static bool reportSaveTime;
	static bool reportTotalLoadTime;
	static bool reportTotalSaveTime;
	static double totalLoadTime;
	static double totalSaveTime;
	static bool trackMemory;
	static bool reportCreation;
	static bool reportManager;
	static bool reportFileBytes;
	static bool countFileBytes;
	static std::uint64_t totalFileReads;
	static std::uint64_t totalFileBytesRead;

public:
	static std::string strip_filter_suffix
	    (const std::string& filename, const int complete=3);
	static std::string default_filter_name
	    (const std::string& filename, const int componentNumber=-1);
	static std::string filter_kind_to_string(std::uint32_t bfKind, bool shortString=true);
	static int vectors_per_filter (const std::uint32_t bfKind);
	static BloomFilter* bloom_filter
	    (const std::string& filename);
	static BloomFilter* bloom_filter
	    (const std::uint32_t bfKind,
	     const std::string& filename, uint32_t kmerSize,
	     std::uint32_t numHashes, std::uint64_t hashSeed1, std::uint64_t hashSeed2,
	     std::uint64_t numBits, std::uint64_t hashModulus=0);
	static BloomFilter* bloom_filter
	    (const BloomFilter* templateBf, const std::string& newFilename="");
	static std::vector<std::pair<std::string,BloomFilter*>> identify_content
	    (std::ifstream& in, const std::string& filename="");
	static double false_positive_rate
		(std::uint32_t numHashes,std::uint64_t numBits,std::uint64_t numItems);

public:
	bool dbgBV               = false;	// some of these are only meaningful
	bool dbgAdd              = false;	// .. if bloom_filter_supportDebug is
	bool dbgContains         = false;	// .. #defined in bloom_filter.cc
	bool dbgAdjustPosList    = false;
	bool dbgRankSelectLookup = false;

public:
	static const std::uint64_t npos = -1;
	static const int unresolved = -1;  // results for lookup()
	static const int absent     =  0;
	static const int present    =  1;
	};


class AllSomeFilter: public BloomFilter
	{
	// nota bene:
	//   numBitVectors is 2
	//   bvAll  is bvs[0]
	//   bvSome is bvs[1]

public:
	AllSomeFilter(const std::string& filename);
	AllSomeFilter(const std::string& filename, uint32_t kmerSize,
	              std::uint32_t numHashes, std::uint64_t hashSeed1, std::uint64_t hashSeed2,
	              std::uint64_t numBits, std::uint64_t hashModulus=0);
	AllSomeFilter(const BloomFilter* templateBf, const std::string& newFilename="");
	virtual ~AllSomeFilter();

	virtual std::string class_identity() const { return "AllsomeFilter"; }
	virtual std::uint32_t kind() const { return bfkind_allsome; }

	virtual void add (const std::string& mer);
	virtual void add (const std::uint64_t* merData);
	virtual bool contains (const std::string& mer) const;
	virtual bool contains (const std::uint64_t* merData) const;
	virtual int lookup (const std::uint64_t pos) const;
	};


class DeterminedFilter: public BloomFilter
	{
	// nota bene:
	//   numBitVectors is 2
	//   bvDetermined is bvs[0]
	//   bvHow        is bvs[1]

public:
	DeterminedFilter(const std::string& filename);
	DeterminedFilter(const std::string& filename, uint32_t kmerSize,
	              std::uint32_t numHashes, std::uint64_t hashSeed1, std::uint64_t hashSeed2,
	              std::uint64_t numBits, std::uint64_t hashModulus=0);
	DeterminedFilter(const BloomFilter* templateBf, const std::string& newFilename="");
	virtual ~DeterminedFilter();

	virtual std::string class_identity() const { return "DeterminedFilter"; }
	virtual std::uint32_t kind() const { return bfkind_determined; }

	virtual int lookup (const std::uint64_t pos) const;
	};

class DeterminedBriefFilter: public DeterminedFilter
	{
	// nota bene:
	//   numBitVectors is 2
	//   bvDetermined is bvs[0]
	//   bvHow        is bvs[1]

	// values for filterInfo; we use filterInfo to distinguish squeezed bit
	// vectors from not-yet-squeezed bit vectors 
	// $$$ can/should I check this at load time?  I only expect to load a
	//     .. not-yet-squeezed bit vector if I'm in the process of building a
	//     .. tree and had to unload/reload a node before it was finished
public:
	enum { squeezed=0, notSqueezed=1 };

public:
	DeterminedBriefFilter(const std::string& filename);
	DeterminedBriefFilter(const std::string& filename, uint32_t kmerSize,
	              std::uint32_t numHashes, std::uint64_t hashSeed1, std::uint64_t hashSeed2,
	              std::uint64_t numBits, std::uint64_t hashModulus=0);
	DeterminedBriefFilter(const BloomFilter* templateBf, const std::string& newFilename="");
	virtual ~DeterminedBriefFilter();

	virtual std::string class_identity() const { return "DeterminedBriefFilter"; }
	virtual std::uint32_t kind() const { return bfkind_determined_brief; }

	virtual int lookup (const std::uint64_t pos) const;

	virtual bool is_position_adjustor  () { return true; }
	virtual void adjust_positions_in_list  (std::vector<std::uint64_t> &kmerPositions,
	                                        std::uint64_t numUnresolved);
	virtual void restore_positions_in_list (std::vector<std::uint64_t> &kmerPositions,
	                                        std::uint64_t numUnresolved);
	};

#endif // bloom_filter_H
