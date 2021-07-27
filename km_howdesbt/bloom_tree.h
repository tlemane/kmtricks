#ifndef bloom_tree_H
#define bloom_tree_H

#include <string>
#include <vector>
#include <iostream>

#include "bloom_filter.h"
#include "query.h"
//#include "kmtricks/repartition.hpp"

class FileManager;

//----------
//
// classes in this module--
//
//----------

struct querystats
	{
	bool			examined;
	bool			passed;
	bool			failed;
	std::uint64_t	numPassed;
	std::uint64_t	numFailed;
	std::uint64_t	numUnresolved;
	std::uint64_t	locallyPassed;
	std::uint64_t	locallyFailed;
	};

enum topofmt
	{
	topofmt_nodeNames = 0,
	topofmt_fileNames,
	topofmt_containers
	};


class BloomTree
	{
public:
	BloomTree(const std::string& name="", const std::string& bfFilename="");
  BloomTree(const std::string& name, const std::string& bfFilename, const std::string& repartFileName, const std::string& winFileName);
	BloomTree(BloomTree* tree);
	virtual ~BloomTree();

	virtual void preload();
	virtual void load();
	virtual void save(bool finished=true);
	virtual void unloadable();

	// virtual void relay_debug_settings();

	virtual void add_child(BloomTree* offspring);
	virtual void disown_children();
	virtual size_t num_children() { return children.size(); }
	virtual bool is_dummy() const { return isDummy; }
	virtual bool is_leaf() const { return isLeaf; }
	virtual BloomTree* child(size_t childNum);
	virtual BloomFilter* real_filter();

	virtual void pre_order (std::vector<BloomTree*>& order);
	virtual void post_order (std::vector<BloomTree*>& order);
	virtual void leaves (std::vector<BloomTree*>& order);

	virtual void print_topology (std::ostream& out, int level=0, int format=topofmt_fileNames) const;
	virtual void construct_union_nodes (std::uint32_t compressor);
	virtual void construct_allsome_nodes (std::uint32_t compressor);
	virtual void construct_determined_nodes (std::uint32_t compressor);
	virtual void construct_determined_brief_nodes (std::uint32_t compressor);
	virtual void construct_intersection_nodes (std::uint32_t compressor);

	virtual void batch_query (std::vector<Query*> queries,
	                          bool completeKmerCounts=false);
private:
	virtual void perform_batch_query (std::uint64_t activeQueries, std::vector<Query*> queries,
	                                  bool completeKmerCounts=false);
	virtual void query_matches_leaves (Query* q);

public:
	virtual int lookup (const std::uint64_t pos) const;
	virtual void enable_query_stats(const std::uint32_t batchSize);
	virtual void clear_query_stats(querystats& stats);
	virtual bool report_query_stats(std::ostream& s,Query* q,bool quietly=true);

public:
	bool isDummy;						// a dummy has no filter; the root might
										// .. be a dummy, to allow for forests
	FileManager* manager;
	std::string name;
	std::string bfFilename;
	std::string futureBfFilename;
	BloomFilter* bf;
	bool isLeaf;
	BloomTree* parent;
	std::vector<BloomTree*> children;	// this will either be empty or have size
										// .. at least 2 (never size 1)
	bool fpRateKnown;
	double fpRate;						// bloom filter false positive rate

	bool nodesShareFiles;				// (only applicable at root)
										// true => tree may contain nodes that
										//         .. share files with each other

public:
	bool reportLoad = false;
	bool reportSave = false;
	static bool inhibitBvSimplify;
	static bool trackMemory;
	static bool reportUnload;
	static int  dbgTraversalCounter;

public:
	std::uint32_t queryStatsLen;
	querystats* queryStats;				// querystats queryStats[queryStatsLen]

public:
	std::uint32_t depth;				// object variables for use by "user"
	std::uint32_t height;				// .. processes
	std::uint32_t subTreeSize;

public:
	bool dbgTraversal           = false;
	bool dbgSortKmerPositions   = false;
	bool dbgKmerPositions       = false;
	// bool dbgKmerPositionsByHash = false;
	// bool dbgLookups             = false;
	// bool dbgInhibitChildUpdate  = false;
	// bool dbgRankSelectLookup    = false;

public:
	std::string repartFileName;
	std::string winFileName;
public:
	static BloomTree* read_topology(const std::string& filename);
	};

#endif // bloom_tree_H
