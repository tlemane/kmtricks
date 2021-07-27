#pragma once

#include <string>
#include <cstdlib>
#include <cstdint>
#include <iostream>

#include "query.h"
#include "commands.h"

class QueryCommandKm: public Command
{
public:
  static constexpr double defaultQueryThreshold = 0.7;  // 70%

public:
  QueryCommandKm(const std::string& name): Command(name) {}
  virtual ~QueryCommandKm();
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
  std::vector<std::string> queryFilenames;
  std::vector<double> queryThresholds;
  std::string matchesFilename;
  std::string repartFileName;
  std::string winFileName;
  double generalQueryThreshold;
  bool sortByKmerCounts;
  bool useFileManager;
  bool checkConsistency;			// only meaningful if useFileManager is false
  // bool collectNodeStats;
  bool reportTime;
  bool completeKmerCounts;

  std::vector<Query*> queries;
};

