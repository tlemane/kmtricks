#ifndef cmd_query_H
#define cmd_query_H

#include <string>
#include <cstdlib>
#include <cstdint>
#include <iostream>

#include "query.h"
#include "commands.h"

class QueryCommand: public Command
	{
public:
	static constexpr double defaultQueryThreshold = 0.7;  // 70%

public:
	QueryCommand(const std::string& name): Command(name) {}
	virtual ~QueryCommand();
	virtual void short_description (std::ostream& s);
	virtual void usage (std::ostream& s, const std::string& message="");
	virtual void debug_help (std::ostream& s);
	virtual void parse (int _argc, char** _argv);
	virtual int execute (void);
	virtual void read_queries (void);
	virtual void sort_matches_by_kmer_counts (void);
	virtual void print_matches(std::ostream& out) const;
	virtual void print_matches_with_kmer_counts(std::ostream& out) const;
	virtual void print_kmer_hit_counts(std::ostream& out) const;

	std::string treeFilename;
    std::string repartFileName;
    std::string winFileName;
	std::vector<std::string> queryFilenames;
	std::vector<double> queryThresholds;
	std::string matchesFilename;
	double generalQueryThreshold;
	bool sortByKmerCounts;
	bool distinctKmers;
	bool useFileManager;
	bool checkConsistency;			// only meaningful if useFileManager is false
	bool justReportKmerCounts;
	bool countAllKmerHits;
	bool reportNodesExamined;
	bool reportTime;
	bool completeKmerCounts;

	std::vector<Query*> queries;
	};

#endif // cmd_query_H
