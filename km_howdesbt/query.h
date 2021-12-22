#ifndef query_H
#define query_H

#include <unordered_set>
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


template<size_t KSIZE>
struct KmerHash
{
  void operator()(const std::string& mer, std::shared_ptr<km::HashWindow> hw, std::shared_ptr<km::Repartition> repart, uint32_t  minim, uint64_t& pos)
  {
    km::Kmer<KSIZE> kmer(mer); km::Kmer<KSIZE> cano = kmer.canonical();
    uint32_t part = repart->get_partition(cano.minimizer(minim).value());
    pos = km::KmerHashers<1>::WinHasher<KSIZE>(part, hw->get_window_size_bits())(cano);
  }
};

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

	virtual void smerize (BloomFilter* bf);

public:
	std::uint32_t batchIx;	// index of this query within a batch
	std::string name;
	std::string seq;		// nucleotide sequence
	double threshold;		// search threshold
	std::vector<std::pair<std::uint64_t,std::size_t>> smerHashes; // the <smers, position> 
										// .. (with smer converted to hash
										// .. values) corresponding to this
										// .. query; the first numUnresolved
										// .. entries are the yet-to-be-resolved
										// .. smers; the resolved smers are
										// .. moved to the tail
	std::vector <size_t> pos_present_smers ; // contains  positions 
										// .. of positive smers at a 
										// .. given node of the tree
										// $$$ think up a destructor
	std::vector<std::unordered_set <size_t>> pos_present_smers_stack ; // for each leave 
										// .. that matches this query,
										// .. store positive positions of kmer set
										// $$$ think up a destructor
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
	std::uint64_t numHashes;            // total size of smerHashes
	
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
