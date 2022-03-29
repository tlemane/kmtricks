#ifndef cmd_query_H
#define cmd_query_H

#include <string>
#include <cstdlib>
#include <cstdint>
#include <iostream>

#include "query.h"
#include "commands.h"


#include <kmtricks/loop_executor.hpp>

class QueryCommand: public Command
	{
public:
	static constexpr double defaultQueryThreshold = 0.7;  // 70%

public:
	QueryCommand(const std::string& name): Command(name) {}
	virtual ~QueryCommand();
	virtual void short_description (std::ostream& s);
	virtual void usage (std::ostream& s, const std::string& message="");
	virtual void parse (int _argc, char** _argv);
	virtual int execute (void);
	virtual void read_queries (void);
	std::vector<bool> get_positive_kmers(const std::string& sequence, 
											const std::unordered_set<std::size_t>& local_presentHashes, 
											const unsigned int& smerSize) const;
	virtual void print_matches_with_kmer_counts_and_spans
   (	std::ostream& out,
		const unsigned int& smerSize
   ) const;
	

	std::string treeFilename;
    std::string repartFileName;
    std::string winFileName;
	std::vector<std::string> queryFilenames;
	std::vector<float> queryThresholds;
	std::string matchesFilename;
	float generalQueryThreshold;
	float threshold_shared_positions;
	bool nodetail;
	bool useFileManager;
	bool checkConsistency;			// only meaningful if useFileManager is false
	bool completeSmerCounts;
	int z; 							// findere strategy

	// needed for findere approach: from smers to hash values when printing results
    std::shared_ptr<km::Repartition> repartitor; 
    std::shared_ptr<km::HashWindow> hash_win; 

	std::vector<Query*> queries;
	};

#endif // cmd_query_H
