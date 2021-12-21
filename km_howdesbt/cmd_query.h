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
	virtual void sort_matches_by_smer_counts (void);
	virtual void print_matches(std::ostream& out) const;
	virtual void print_matches_with_smer_counts(std::ostream& out) const;
	virtual std::vector<bool> positiveKmers(const std::string& sequence, const std::unordered_set<std::uint64_t>& local_presentHashes, const unsigned int& smerSize, const unsigned int& z) const ;

	std::string treeFilename;
    std::string repartFileName;
    std::string winFileName;
	std::vector<std::string> queryFilenames;
	std::vector<double> queryThresholds;
	std::string matchesFilename;
	double generalQueryThreshold;
	bool sortBySmerCounts;
	bool useFileManager;
	bool checkConsistency;			// only meaningful if useFileManager is false
	bool completeSmerCounts;
	int z; 							// findere strategy

	std::vector<Query*> queries;
	};

#endif // cmd_query_H
