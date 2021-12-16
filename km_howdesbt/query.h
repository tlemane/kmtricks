#ifndef query_H
#define query_H

#include <string>
#include <vector>
#include <iostream>

#include "bloom_filter.h"

#define KMTRICKS_PUBLIC
#define WITH_XXHASH
#include <kmtricks/kmer.hpp>
#include <kmtricks/kmer_hash.hpp>
#include <kmtricks/hash.hpp>
#include <kmtricks/repartition.hpp>
#include <kmtricks/loop_executor.hpp>

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
    Query(const querydata& qd, double threshold, std::shared_ptr<km::Repartition> rep, std::shared_ptr<km::HashWindow> hash_win, uint32_t minimsize);
    virtual ~Query();

	virtual void smerize (BloomFilter* bf, bool distinct=false, bool populateSmers=false);
	virtual void sort_smer_positions ();
	virtual void dump_smer_positions (std::uint64_t numUnresolved=-1);
	virtual std::uint64_t smer_positions_hash (std::uint64_t numUnresolved=-1);

public:
	std::uint32_t batchIx;	// index of this query within a batch
	std::string name;
	std::string seq;		// nucleotide sequence
	double threshold;		// search threshold
	std::vector<std::uint64_t> smerHashes; // the smers (converted to hash
										// .. values) corresponding to this
										// .. query; the first numUnresolved
										// .. entries are the yet-to-be-resolved
										// .. smers; the resolved smers are
										// .. moved to the tail
	std::vector<std::string> smers;		// the smers; this is only populated
										// .. in special instances (e.g. for
										// .. cmd_query_bf), and in those
										// .. cases care should be taken to
										// .. assure that the smerHashes[ix]
										// .. corresponds to smers[ix] for each
										// .. ix

	std::uint64_t numHashes;			// total size of smerHashes
	std::uint64_t neededToPass;			// number of smers required, to judge
										// .. the query as a "pass"
	std::uint64_t neededToFail;			// number of smers required, to judge
										// .. the query as a "fail"
	std::uint64_t numUnresolved;		// number of smers not yet known to be
										// .. present or absent in all leaves
										// .. of the current subtree
	std::uint64_t numPassed;			// number of smers known to be present
										// .. in all leaves of the subtree
	std::uint64_t numFailed;			// number of smers known to be absent
										// .. in all leaves of the subtree
	bool adjustsmerCounts;				// true  => populate matchesAdjusted[]
										// false => don't
    std::vector<std::string> matches;	// names of leaves that match this query
    std::vector<std::uint64_t> matchesNumPassed;  // numPassed corresponding to
										// .. each match; only valid if the
										// .. search reached the leaf without
										// .. having been pruned
    std::vector<std::uint64_t> matchesAdjustedHits;  // adjusted value of
										// .. numPassed (corresponding to each
										// .. match), to .. account for
										// .. estimated bloom filter false
										// .. positives; only valid if the
										// .. search reached the leaf without
										// .. having been pruned, and only if
										// .. adjustsmerCounts is true

										// stacks to maintain numUnresolved,
										// .. numPassed, and numFailed as we
										// .. move up and down the tree
    std::vector<std::uint64_t> numUnresolvedStack;
    std::vector<std::uint64_t> numPassedStack;
    std::vector<std::uint64_t> numFailedStack;

    std::shared_ptr<km::Repartition> m_repartitor {nullptr};
    std::shared_ptr<km::HashWindow> m_hash_win {nullptr};
    uint64_t m_m_size;
    uint64_t m_minim_size;

public:
	static void read_query_file (std::istream& in, const std::string& filename,
	                             double threshold,
	                             std::vector<Query*>& queries,
                                 std::string& repartFileName, std::string& winFileName);
	};

#endif // query_H
