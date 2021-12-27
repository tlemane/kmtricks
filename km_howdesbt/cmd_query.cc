// cmd_query.cc-- query a sequence bloom tree

#include <string>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <tuple>

#include "utilities.h"
#include "bit_vector.h"
#include "bloom_filter.h"
#include "bloom_tree.h"
#include "file_manager.h"
#include "query.h"

#include "support.h"
#include "commands.h"
#include "cmd_query.h"

using std::string;
using std::vector;
using std::pair;
using std::tuple;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
#define u32 std::uint32_t
#define u64 std::uint64_t


void QueryCommand::short_description
   (std::ostream& s)
	{
	s << commandName << "-- query a sequence bloom tree" << endl;
	}

void QueryCommand::usage
   (std::ostream& s,
	const string& message)
	{
	if (!message.empty())
		{
		s << message << endl;
		s << endl;
		}

//$$$ add an option to limit the number of bits used in each BF
//$$$ .. that's to let us experiment with different reductions of BF fraction
//$$$ .. without having to generate every populated filter size; implementation
//$$$ .. would just act as a filter on the hashed position list for each query
	short_description(s);
	s << "usage: " << commandName << " [<queryfilename>[=<F>]] [options]" << endl;
	//    123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	s << "  --tree=<filename>    name of the tree toplogy file" << endl;
	s << "  <queryfilename>      (cumulative) name of a query file; this is either a" << endl;
	s << "                       fasta file or a file with one nucleotide sequence per" << endl;
	s << "                       line; if no query files are provided, queries are read" << endl;
	s << "                       from stdin" << endl;
	s << "  <queryfilename>=<F>  query file with associated threshold; <F> has the same" << endl;
	s << "                       meaning as in --threshold=<F> but applies only to this" << endl;
	s << "                       query file" << endl;
	s << "  --repart=<F>         minimizers repartition (from kmtricks)" << endl;
    s << "  --win=<F>            hash window (from kmtricks)" << endl;
    s << "  --threshold=<F>      fraction of query kmers that must be present in a leaf" << endl;
	s << "                       to be considered a match; this must be between 0 and 1;" << endl;
	s << "                       this only applies to query files for which <F> is not" << endl;
	s << "                       otherwise specified (by <queryfilename>=<F>)" << endl;
	s << "                       (default is " << defaultQueryThreshold << ")" << endl;
	s << "  --threshold-shared-positions=<F> Prints a query result if its ratio" << endl;
    s << "                       of positions covered by at least a shared kmer is" << endl;
    s << "                       higher or equal to this threshold. This happens" << endl;
    s << "                       after the threshold applied on the" << endl;
    s << "                       number of shared kmers. This option enables to" << endl;
    s << "                       save query results where, say, 60 of kmers are" << endl;
    s << "                       shared but 95% of positions are covered by a " << endl;
    s << "                       shared kmer. In this case with this value set to 90, " << endl;
    s << "                       this result is printed." << endl;
	s << "                       (default is " << defaultQueryThreshold << ")" << endl;
	s << "  --no-detail          Do not print the position of shared kmers in output." << endl;
	s << "	--z=<F>				 If z is bigger than 0, apply the findere strategy." << endl;
	s << " 						 In such case, a k-mer is considered as present" <<endl;
	s << " 						 if all its s-mers are presents," << endl;
	s << " 						 with k = s+z, and s being the size of the words indexed in"<< endl;
	s << " 						 bloom filters. Hence, with z=0 (defaut value), no fidere " <<endl;
	s << " 						 approach is applied, and words indexed in the blomm filters" <<endl;
	s << " 						 are queried" << endl; 
	s << "  --consistencycheck   before searching, check that bloom filter properties are" << endl;
	s << "                       consistent across the tree" << endl;
	s << "                       (not needed with --usemanager)" << endl;
	s << "  --time               report wall time and node i/o time" << endl;
	s << "  --out=<filename>     file for query results; if this is not provided, results" << endl;
	s << "                       are written to stdout" << endl;

	}
void QueryCommand::parse
   (int		_argc,
	char**	_argv)
	{
	int		argc;
	char**	argv;

	// defaults

	generalQueryThreshold   	= -1.0;		// (unassigned threshold)
	nodetail        			= false;
	threshold_shared_positions 	= defaultQueryThreshold;
	checkConsistency        	= false;
	z							= 0;


	// skip command name

	argv = _argv+1;  argc = _argc - 1;
	if (argc <= 0) chastise ();

	//////////
	// scan arguments
	//////////

	for (int argIx=0 ; argIx<argc ; argIx++)
		{
		string arg = argv[argIx];
		string argVal;
		if (arg.empty()) continue;

		string::size_type argValIx = arg.find('=');
		if (argValIx == string::npos) argVal = "";
		else argVal = arg.substr(argValIx+1);

		// --help, etc.
		if ((arg == "--help")
		 || (arg == "-help")
		 || (arg == "--h")
		 || (arg == "-h")
		 || (arg == "?")
		 || (arg == "-?")
		 || (arg == "--?"))
			{ usage (cerr);  std::exit (EXIT_SUCCESS); }


		// --tree=<filename>, etc.

		if ((is_prefix_of (arg, "--tree="))
		 ||	(is_prefix_of (arg, "--intree="))
		 ||	(is_prefix_of (arg, "--topology=")))
			{ treeFilename = argVal;  continue; }

		// (unadvertised) --query=<filename>
		//             or --query=<filename>=<F> or --query=<filename>:<F>

		if (is_prefix_of (arg, "--query="))
			{
			string::size_type threshIx = argVal.find('=');
			if (threshIx == string::npos) threshIx = argVal.find(':');

			if (threshIx == string::npos)
				{
				queryFilenames.emplace_back(strip_blank_ends(argVal));
				queryThresholds.emplace_back(-1.0);  // (unassigned threshold)
				}
			else
				{
				double thisQueryThreshold = string_to_probability(arg.substr(threshIx+1));
				queryFilenames.emplace_back(strip_blank_ends(argVal));
				queryThresholds.emplace_back(thisQueryThreshold);
				}
			continue;
			}

        if (is_prefix_of(arg, "--repart="))
        {
			repartFileName = argVal;
			continue;
        }

        if (is_prefix_of(arg, "--win="))
        {
			winFileName = argVal;
			continue;
        }

		// --z=<F>
		if (is_prefix_of(arg, "--z="))
        {
			z = std::stod (argVal);
			continue;
        }

		// --threshold=<F>

		if ((is_prefix_of (arg, "--threshold="))
		 ||	(is_prefix_of (arg, "--query-threshold="))
		 ||	(is_prefix_of (arg, "--theta="))
		 ||	(is_prefix_of (arg, "--specificity=")))
			{
			if (generalQueryThreshold >= 0.0)
				{
				cerr << "warning: --threshold=<F> used more that once; only final setting will apply" << endl;
				cerr << "(to use different thresholds for different files, use <queryfilename>=<F> form)" << endl;
				}
			generalQueryThreshold = string_to_probability(argVal);
			continue;
			}

		// --threshold-shared-positions=<F>

		if ((is_prefix_of (arg, "--threshold-shared-positions="))
		 ||	(is_prefix_of (arg, "--threshold_shared_positions=")))
			{
				threshold_shared_positions = string_to_probability(argVal);
				continue;
			}

		// --no-detail

		if (arg == "--no-detail")
			{ nodetail = true;  continue; }




		// --consistencycheck, (unadvertised) --noconsistency

		if (arg == "--consistencycheck")
			{ checkConsistency = true;  continue; }

		if ((arg == "--noconsistency")
		 || (arg == "--noconsistencycheck"))
			{ checkConsistency = false;  continue; }



		
	
		// --out=<filename>, etc.

		if ((is_prefix_of (arg, "--out="))
		 ||	(is_prefix_of (arg, "--output="))
		 ||	(is_prefix_of (arg, "--matches="))
		 ||	(is_prefix_of (arg, "--results=")))
			{ matchesFilename = argVal;  continue; }

		// (unadvertised) debug options

		if (arg == "--debug")
			{ debug.insert ("debug");  continue; }

		if (is_prefix_of (arg, "--debug="))
			{
		    for (const auto& field : parse_comma_list(argVal))
				debug.insert(to_lower(field));
			continue;
			}

		// unrecognized --option

		if (is_prefix_of (arg, "--"))
			chastise ("unrecognized option: \"" + arg + "\"");

		// <queryfilename>=<F> or <queryfilename>:<F>

		string::size_type threshIx = argValIx;
		if (threshIx == string::npos) threshIx = arg.find(':');

		if (threshIx != string::npos)
			{
			double thisQueryThreshold = string_to_probability(arg.substr(threshIx+1));
			queryFilenames.emplace_back(strip_blank_ends(arg.substr(0,threshIx)));
			queryThresholds.emplace_back(thisQueryThreshold);
			continue;
			}

		// <queryfilename>

		queryFilenames.emplace_back(strip_blank_ends(arg));
		queryThresholds.emplace_back(-1.0);  // (unassigned threshold)
		}

	// sanity checks

	if (treeFilename.empty())
		chastise ("you have to provide a tree topology file");



	repartitor = std::make_shared<km::Repartition>(repartFileName, "");
	hash_win = std::make_shared<km::HashWindow>(winFileName);

	completeSmerCounts = true; // $$$ remove this ?

	// assign threshold to any unassigned queries

	if (generalQueryThreshold < 0.0)
		generalQueryThreshold = defaultQueryThreshold;

	int numQueryFiles = queryFilenames.size();
	for (int queryIx=0 ; queryIx<numQueryFiles ; queryIx++)
		{
		if (queryThresholds[queryIx] < 0)
			queryThresholds[queryIx] = generalQueryThreshold;
		}

	return;
	}

QueryCommand::~QueryCommand()
	{
	for (const auto& q : queries)
		delete q;
	}

int QueryCommand::execute()
	{


	// read the tree

	BloomTree* root = BloomTree::read_topology(treeFilename);

	useFileManager = root->nodesShareFiles;

	vector<BloomTree*> order;





	// set up the file manager

	FileManager* manager = nullptr;
	if (useFileManager)
		{
		manager = new FileManager(root,/*validateConsistency*/false);
		}

	// if we're not using a file manager, we may still want to do a consistency
	// check before we start the search (we'd rather not run for a long time
	// and *then* report the problem)

	else if (checkConsistency)
		{
		BloomFilter* modelBf = nullptr;

		if (order.size() == 0)
			root->post_order(order);
		for (const auto& node : order)
			{
			node->preload();

			if (modelBf == nullptr)
				modelBf = node->bf;
			else
				node->bf->is_consistent_with (modelBf, /*beFatal*/ true);
			}
		}
	

	// read the queries

	read_queries ();


	// perform the query

	root->batch_query(queries,completeSmerCounts);


	// get the smer size
	// TODO dirty, how to easily get the s value?
	unsigned int smerSize = 0;
	BloomFilter* modelBf = nullptr;

	if (order.size() == 0)
		root->post_order(order);
	for (const auto& node : order)
		{
		node->preload();
		smerSize = node->bf->smerSize;
		break;
		}			
		

	


	std::streambuf * buf;
	std::ofstream of;

	if(!matchesFilename.empty()) {
    	of.open(matchesFilename);
    	buf = of.rdbuf();
	} 
	else {
    	buf = std::cout.rdbuf();
	}
	std::ostream out(buf);


	
	
	print_matches_with_kmer_counts_and_spans (out, smerSize);
	
			
	
		

//$$$ where do we delete the tree?  looks like a memory leak

	FileManager::close_file();	// make sure the last bloom filter file we
								// .. opened for read gets closed

	if (manager != nullptr)
		delete manager;





	return EXIT_SUCCESS;
	}

//----------
//
// read_queries--
//	Read the query file(s), populating the queries list.
//
// $$$ this should warn the user if the queries have any name used more than
//     once for sequences that aren't the same
//----------

void QueryCommand::read_queries()
	{
	// if no query files are provided, read from stdin

	if (queryFilenames.empty())
		Query::read_query_file (cin, /*filename*/ "", generalQueryThreshold, queries, repartFileName, winFileName);

	// otherwise, read each query file

	else
		{
		int numQueryFiles = queryFilenames.size();
		for (int queryIx=0 ; queryIx<numQueryFiles ; queryIx++)
			{
			string filename = queryFilenames[queryIx];
			std::ifstream in (filename);
			if (not in)
				fatal ("error: failed to open \"" + filename + "\"");
			Query::read_query_file (in, filename, queryThresholds[queryIx], queries, repartFileName, winFileName);
			in.close();
			}
		}

	}





/**
 * @brief From hash values associated to smers to a vector of Positive kmers. Code adapted from findere https://github.com/lrobidou/findere/, by Lucas Robidou
 * @param sequence: input sequence.
 * @param local_presentHashes: s-mer hash values positive for this sequence
 * @param smerSize size of indexed s-mers
 * @param z determines k, the  size of queried k-mers (k=s+z)
 * @return The result of findere's query on the sequence.
 */
std::vector<bool> QueryCommand::get_positive_kmers(const std::string& sequence, 
											const std::unordered_set<std::size_t>& local_pos_present_smers, 
											const unsigned int& smerSize) const {
	unsigned long long size = sequence.size();
    const unsigned int kmerSize = smerSize + z; 
    std::vector<bool> response(size - kmerSize + 1, false);
	if (z == 0){
		for (unsigned long long j = 0; j < size - smerSize + 1; j++)
			response[j] = local_pos_present_smers.count(j) > 0;
		return response;
	}
    unsigned long long stretchLength = 0;  // number of consecutive positives kmers
    unsigned long long j = 0;              // index of the query vector
    bool extending_stretch = true;

    while (j < size - smerSize + 1) 
	{
		if (local_pos_present_smers.count(j) > 0) {
            if (extending_stretch) {
                stretchLength++;
                j++;
            } else {
                extending_stretch = true;
                j = j - z;
            }
        } else {
            if (stretchLength >= z) {
                for (unsigned long long t = j - stretchLength; t < j - z; t++) response[t] = true;
            }
            stretchLength = 0;
            extending_stretch = false;
            j = j + z + 1;
        }
    }
    // Last values:
    if (stretchLength >= z) {
        for (unsigned long long t = size - smerSize + 1 - stretchLength; t < size - kmerSize + 1; t++) response[t] = true;
    }

    return response;
}


/**
 * @brief Computes a boolean vector of positions covered by at least a shared kmer.
 * @param shared_kmers The result of the query.
 * @param kmerSize the value of k.
 * @return a boolean vector of positions covered by at least a shared kmer.
 */
std::vector<bool> get_positions_covered(std::vector<bool> shared_kmers, const unsigned int kmerSize) {
	// assumes that shared_kmers.size() > 0
	std::vector<bool> response(shared_kmers.size() + kmerSize - 1, false);
    long long last_covered_position = -1;
    long long pos = 0;
    for (auto shared : shared_kmers) {
        if (shared) {
            if (last_covered_position < pos) {
				for (int i = pos; i < pos + kmerSize; i++)
					response[i] = true;
            } else {
				for (int i = last_covered_position + 1; i < pos + kmerSize; i++)
					response[i] = true;
            }
            last_covered_position = pos + kmerSize - 1;
        }
        pos++;
    }
    return response;
}


//----------
//
// print_matches_with_kmer_counts_and_spans--
//
//----------
void QueryCommand::print_matches_with_kmer_counts_and_spans
   (	std::ostream& out,
		const unsigned int& smerSize
   ) const
	{
	std::ios::fmtflags saveOutFlags(out.flags());

	
	out << "# FORMAT:" << endl;
	out << "# * [query name]" <<endl;
	out << "# For each target, 3 or 4 fields:" <<endl;
	out << "#   [taget name]" <<endl;
	out << "#   (unless --no-detail option) string in {+-} showing positions covered (+) by at least a shared kmer, else (-)" <<endl;
	out << "#   Ratio of kmers of the query shared with the target"  <<endl;
	out << "#   Ratio of positions of the query covered by at least a kmer shared with the target"  <<endl;


	for (auto& q : queries)
		{
		out << "* [" << q->name << "] "<< endl;
		std::string   seq = q->seq;
		// For each query, we store its answers in a vector of tuples:
		// .. <float, string, float, string>, being: 
		// .. 0 <ratio position covered by a shared kmer (key of sorting), 
		// .. 1 name of the target reference, 
		// .. 2 ratio nb shared kmers, 
		// .. 3 +-. string indicating the position of the shared kmers in the copy>
		std::vector<std::tuple<float, string, float, string>> res_matches;
		int matchIx = 0;
		for (auto& name : q->matches)
			{
			u64 numPassed = q->matchesNumPassed[matchIx];
			std::vector<bool> positive_kmers = get_positive_kmers(	seq, 
																	q->pos_present_smers_stack[matchIx], 
																	smerSize );

			std::vector<bool> positions_covered = get_positions_covered(positive_kmers, smerSize + z);
			unsigned long long nb_positions_covered = std::count(positions_covered.begin(), positions_covered.end(), true);
			
			std::string pmres = "";
			bool pm_shows_ratio_kmers = false; // $$$ remove when chosen.
			if (not nodetail)
			{
				if (pm_shows_ratio_kmers) // prints position starting a shared kmer
				{
				for (auto is_present: positive_kmers)
					{
					if (is_present) pmres.append("+");
					else pmres.append("-");
					}
				for (int i = 0; i < smerSize + z - 1; i++) pmres.append("."); // last kmer 
				}
				else // prints position covered by at least a shared kmer
				{
				for (auto is_present: positions_covered)
					{
					if (is_present) pmres.append("+");
					else pmres.append("-");
					}
				}
				pmres.append(" "); // add a space to simplify the printing when nodetail is required
			}		

			float positive_kmer_ratio = std::count(positive_kmers.begin(), positive_kmers.end(), true)/float(positive_kmers.size());
			float positive_covered_pos_ratio = nb_positions_covered/float(positive_kmers.size() + smerSize + z -1 );
			

			// kmer number <= smer number. Hence we recheck that the kmer threshold does not get below the 
			// required threshold.
			if (positive_kmer_ratio >= q->threshold or positive_covered_pos_ratio >= threshold_shared_positions)
				res_matches.push_back(std::make_tuple(positive_covered_pos_ratio, name, positive_kmer_ratio, pmres));

			matchIx++;
			}
		// sort and print results:
		sort(res_matches.begin(), res_matches.end(), std::greater <>());
		for (auto match: res_matches)
			{
			out << "[" << std::get<1>(match) << "] " << std::get<3>(match);
			out << std::setprecision (2) << std::fixed << std::get<2>(match) << " ";
			out << std::setprecision (2) << std::fixed << std::get<0>(match) << endl;
			}
		}



	out.flags(saveOutFlags);
	}
