#ifndef query_H
#define query_H

#include <string>
#include <vector>
#include <iostream>

#include "bloom_filter.h"
//#include "query.h"
#include "sabuhash.h"

//#define _KM_LIB_INCLUDE_
//#include "kmtricks/kmlib.hpp"
#include <kmtricks/config.hpp>
#include <kmtricks/loop_executor.hpp>
#include <kmtricks/hash.hpp>
#include <kmtricks/kmer.hpp>
#define WITH_XXHASH
#include <kmtricks/kmer_hash.hpp>
#include <kmtricks/repartition.hpp>
//----------
//
// classes in this module--
//
//----------

struct querydata
	{
	std::uint32_t batchIx;	// index of this query within a batch
	std::string   name;
	std::string   seq;		// nucleotide sequence
	};


class Query
	{
public:
	Query(const querydata& qd, double threshold);
	Query(const querydata& qd, double threshold, std::shared_ptr<km::Repartition> rep, std::shared_ptr<km::HashWindow> hashwin, uint32_t minimsize);
	virtual ~Query();

	virtual void kmerize (BloomFilter* bf, bool distinct=false);

public:
	std::uint32_t batchIx;	// index of this query within a batch
	std::string name;
	std::string seq;		// nucleotide sequence
	std::uint64_t seq_length; 
	double threshold;		// search threshold
	std::vector<std::uint64_t> kmerPositions; // the kmers (converted to hash
										// .. values) corresponding to this
										// .. query; the first numUnresolved
										// .. entries are the yet-to-be-resolved
										// .. kmers; the resolved kmers are
										// .. moved to the tail
	
	std::vector<std::uint64_t> kmerized2endpos; // ending position of each queried 
										// .. kmer stored in kmerPositions. 
										// .. Eg. if first kmer stored in kmerPositions
										// ends position 42, and the second ends position 137, 
										// then kmerized2endpos contains 42 and 137.

									
	std::vector<std::uint64_t> endingPositionSharedKmer; // Ending positions of a shared kmers for a target.


	std::uint64_t numPositions;			// total size of kmerPositions
	std::uint64_t neededToPass;			// number of kmers required, to judge
										// .. the query as a "pass"
	std::uint64_t neededToFail;			// number of kmers required, to judge
										// .. the query as a "fail"
	std::uint64_t numUnresolved;		// number of kmers not yet known to be
										// .. present or absent in all leaves
										// .. of the current subtree
	std::uint64_t numPassed;			// number of kmers known to be present
										// .. in all leaves of the subtree
	std::uint64_t numFailed;			// number of kmers known to be absent
										// .. in all leaves of the subtree
	std::uint64_t nodesExamined;		// number of nodes that were "examined"
										// by this query
    std::vector<std::string> matches;	// names of leaves that match this query
    std::vector<std::uint64_t> matchesNumPassed;  // numPassed corresponding to
										// .. each match; only valid if the
										// .. search reached the leaf without
										// .. having been pruned
	
    std::vector<std::uint64_t> matchesCoveredPos;  // nb of positions covered by
										// a shared kmer, corresponding to
										// .. each match; only valid if the
										// .. search reached the leaf without
										// .. having been pruned
										
    std::vector<std::uint64_t> numUnresolvedStack;
    std::vector<std::uint64_t> numPassedStack;
    std::vector<std::uint64_t> numFailedStack;
    

  std::shared_ptr<km::Repartition> _repartitor {nullptr};
	std::shared_ptr<km::HashWindow> _hash_win {nullptr};
  SabuHash *_h;
  uint64_t _msize;
  uint32_t _minimsize;


public:
	static void read_query_file_km (std::istream& in, const std::string& filename, double threshold, std::vector<Query*>& queries, std::string& repartFileName, std::string& winFileName);
	};

#endif // query_H
