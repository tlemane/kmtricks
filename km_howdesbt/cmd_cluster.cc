// cmd_cluster.cc-- determine a tree topology by clustering bloom filters

#include <string>
#include <cstdlib>
#include <cstdint>
#include <limits>
#include <iostream>
#include <queue>

#include "utilities.h"
#include "bit_utilities.h"
#include "bloom_filter.h"
#include "file_manager.h"

#include "support.h"
#include "commands.h"
#include "cmd_build_sbt.h"
#include "cmd_cluster.h"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
#define u32 std::uint32_t
#define u64 std::uint64_t

void ClusterCommand::short_description
   (std::ostream& s)
	{
	s << commandName << "-- determine a tree topology by clustering bloom filters" << endl;
	}

void ClusterCommand::usage
   (std::ostream& s,
	const string& message)
	{
	if (!message.empty())
		{
		s << message << endl;
		s << endl;
		}

	short_description(s);
	s << "usage: " << commandName << " [options]" << endl;
	//    123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	s << "  --list=<filename> file containing a list of bloom filters to cluster; only" << endl;
	s << "                    filters with uncompressed bit vectors are allowed" << endl;
	s << "  <filename>        same as --list=<filename>" << endl;
	s << "  --out=<filename>  name for tree toplogy file" << endl;
	s << "                    (by default this is derived from the list filename)" << endl;
	s << "  --tree=<filename> same as --out=<filename>" << endl;
	s << "  --nodename=<template> filename template for internal tree nodes" << endl;
	s << "                    this must contain the substring {number}" << endl;
	s << "                    (by default this is derived from the list filename)" << endl;
	s << "  <start>..<end>    interval of bits to use from each filter; the clustering" << endl;
	s << "                    algorithm only considers this subset of each filter's bits" << endl;
	s << "                    (by default we use the first " << defaultEndPosition << " bits)" << endl;
	s << "  --bits=<N>        number of bits to use from each filter; same as 0..<N>" << endl;
	s << "  --cull            remove nodes from the binary tree; remove those for which" << endl;
	s << "                    saturation of determined is more than 2 standard deviations" << endl;
	s << "                    below the mean" << endl;
	s << "                    (this is the default)" << endl;
	s << "  --cull=<Z>sd      remove nodes for which saturation of determined is more" << endl;
	s << "                    than <Z> standard deviations below the mean" << endl;
	s << "  --cull=<S>        remove nodes for which saturation of determined is less" << endl;
	s << "                    than <S>; e.g. <S> can be \"0.20\" or \"20%\"" << endl;
	s << "  --keepallnodes    keep all nodes of the binary tree" << endl;
	s << "  --nocull          (same as --keepallnodes)" << endl;
	s << "  --nobuild         perform the clustering but don't build the tree's nodes" << endl;
	s << "                    (this is the default)" << endl;
	s << "  --build           perform clustering, then build the uncompressed nodes" << endl;
	}

void ClusterCommand::debug_help
   (std::ostream& s)
	{
	s << "--debug= options" << endl;
	s << "  trackmemory" << endl;
	s << "  bvcreation" << endl;
	s << "  interval" << endl;
	s << "  offsets" << endl;
	s << "  console" << endl;
	s << "  bits" << endl;
	s << "  distances" << endl;
	s << "  queue" << endl;
	s << "  mergings" << endl;
	s << "  numbers" << endl;
	s << "  cull" << endl;
	s << "  detratio" << endl;
	s << "  detratiodistrib" << endl;
	}

void ClusterCommand::parse
   (int		_argc,
	char**	_argv)
	{
	int		argc;
	char**	argv;

	// defaults

	startPosition          = 0;
	endPosition            = defaultEndPosition;
	cullNodes              = true;
	deriveCullingThreshold = true;
	cullingThresholdSD     = defaultCullingThresholdSD;
	cullingThreshold       = std::numeric_limits<double>::quiet_NaN();
	renumberNodes          = true;
	inhibitBuild           = true;

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

		if ((arg == "--help=debug")
		 || (arg == "--help:debug")
		 || (arg == "?debug"))
			{ debug_help(cerr);  std::exit (EXIT_SUCCESS); }

		// --list=<filename>

		if (is_prefix_of (arg, "--list="))
			{
			if (not listFilename.empty())
				chastise ("unrecognized option: \"" + arg + "\""
				          "\nbloom filters list file was already given as \"" + listFilename + "\"");
			listFilename = argVal;
			continue;
			}

		// --out=<filename>, --tree=<filename>, etc.

		if ((is_prefix_of (arg, "--out="))
		 || (is_prefix_of (arg, "--output="))
		 || (is_prefix_of (arg, "--tree="))
		 || (is_prefix_of (arg, "--outtree="))
		 || (is_prefix_of (arg, "--topology=")))
			{ treeFilename = argVal;  continue; }

		// --nodename=<template>
		// (and, for backward compatibility, --node=<template>)

		if ((is_prefix_of (arg, "--nodename="))
		 || (is_prefix_of (arg, "--nodenames="))
		 || (is_prefix_of (arg, "--node="))
		 || (is_prefix_of (arg, "--nodes=")))
			{
			nodeTemplate = argVal; 

			std::size_t fieldIx = nodeTemplate.find ("{number}");
			if (fieldIx == string::npos)
				fieldIx = nodeTemplate.find ("{number:1}");
			if (fieldIx == string::npos)
				fieldIx = nodeTemplate.find ("{number:0}");

			// for backward compatibility, we allow {node} instead of {number}
			if (fieldIx == string::npos)
				{
				string field = "{node}";
				fieldIx = nodeTemplate.find (field);
				if (fieldIx != string::npos)
					nodeTemplate.replace(fieldIx,field.length(),"{number}");
				}
			if (fieldIx == string::npos)
				{
				string field = "{node:1}";
				fieldIx = nodeTemplate.find (field);
				if (fieldIx != string::npos)
					nodeTemplate.replace(fieldIx,field.length(),"{number:1}");
				}
			if (fieldIx == string::npos)
				{
				string field = "{node:1}";
				fieldIx = nodeTemplate.find (field);
				if (fieldIx != string::npos)
					nodeTemplate.replace(fieldIx,field.length(),"{number:0}");
				}
			if (fieldIx == string::npos)
				chastise ("--node is required to contain the substring \"{number}\", or a variant of it");

			if (not is_suffix_of(nodeTemplate,".bf"))
				nodeTemplate = nodeTemplate + ".bf";
			continue;
			}

		// --bits=<N>

		if ((is_prefix_of (arg, "--bits="))
		 || (is_prefix_of (arg, "B="))
		 || (is_prefix_of (arg, "--B=")))
			{
			startPosition = 0;
			endPosition   = string_to_unitized_u64(argVal);
			continue;
			}

		// --nocull, --cull
		// --nowinnow, --winnow (unadvertised; for backward compatibility)

		if ((arg == "--nocull")
		 || (arg == "--noculling")
		 || (arg == "--dontcull")
		 || (arg == "--keepallnodes")
		 || (arg == "--nowinnow")
		 || (arg == "--nowinnowing")
		 || (arg == "--dontwinnow"))
			{
			cullNodes              = false;
			deriveCullingThreshold = false;
			cullingThresholdSD     = std::numeric_limits<double>::quiet_NaN();
			cullingThreshold       = std::numeric_limits<double>::quiet_NaN();
			continue;
			}

		if ((arg == "--cull")
		 || (arg == "--culling")
		 || (arg == "--winnow")
		 || (arg == "--winnowing"))
			{
			cullNodes              = true;
			deriveCullingThreshold = true;
			cullingThresholdSD     = defaultCullingThresholdSD;
			cullingThreshold       = std::numeric_limits<double>::quiet_NaN();
			continue;
			}

		// --cull=<Z>sd

		if ((is_suffix_of (arg, "sd"))
		 && ((is_prefix_of (arg, "--cull="))
		  || (is_prefix_of (arg, "--culling="))))
			{
			cullNodes              = true;
			deriveCullingThreshold = true;
			cullingThresholdSD     = string_to_double(strip_suffix(argVal,"sd"));
			cullingThreshold       = std::numeric_limits<double>::quiet_NaN();
			continue;
			}

		// --cull=<S>
		// --winnow=<S> (unadvertised; for backward compatibility)

		if ((is_prefix_of (arg, "--cull="))
		 || (is_prefix_of (arg, "--culling="))
		 || (is_prefix_of (arg, "--winnow="))
		 || (is_prefix_of (arg, "--winnowing=")))
			{
			cullNodes              = true;
			deriveCullingThreshold = false;
			cullingThresholdSD     = std::numeric_limits<double>::quiet_NaN();
			cullingThreshold       = string_to_probability(argVal);
			continue;
			}

		// --norenumber (unadvertised)

		if (arg == "--norenumber")
			{ renumberNodes = false;  continue; }

		// --nobuild, --build

		if ((arg == "--nobuild")
		 || (arg == "--dontbuild"))
			{ inhibitBuild = true;  continue; }

		if (arg == "--build")
			{ inhibitBuild = false;  continue; }

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

		// <start>..<end>

		std::size_t separatorIx = arg.find ("..");
		if (separatorIx != string::npos)
			{
			startPosition = string_to_unitized_u64(arg.substr (0, separatorIx));
			endPosition   = string_to_unitized_u64(arg.substr (separatorIx+2));
			if (endPosition <= startPosition)
				chastise ("bad interval: " + arg + " (end <= start)");
			continue;
			}

		// <filename>

		if (not listFilename.empty())
			chastise ("unrecognized option: \"" + arg + "\""
			          "\nbloom filters list file was already given as \"" + listFilename + "\"");
		listFilename = arg;
		}

	// sanity checks

	if (startPosition % 8 != 0)
		chastise ("the bit interval's start (" + std::to_string(startPosition)
		        + ") has to be a multiple of 8");

	if (listFilename.empty())
		chastise ("you have to provide a file, listing the bloom filters for the tree");

	if (treeFilename.empty())
		{
		string::size_type dotIx = listFilename.find_last_of(".");
		if (dotIx == string::npos)
			treeFilename = listFilename + ".sbt";
		else
			treeFilename = listFilename.substr(0,dotIx) + ".sbt";
		}

	if (nodeTemplate.empty())
		{
		string::size_type dotIx = listFilename.find_last_of(".");
		if (dotIx == string::npos)
			nodeTemplate = listFilename + "{number}.bf";
		else
			nodeTemplate = listFilename.substr(0,dotIx) + "{number}.bf";
		}

	return;
	}


ClusterCommand::~ClusterCommand()
	{
	for (const auto& bf : leafVectors)
		delete bf;
	if (treeRoot != nullptr) delete treeRoot;
	}


int ClusterCommand::execute()
	{
	if (contains(debug,"trackmemory"))
		{
		trackMemory              = true;
		FileManager::trackMemory = true;
		BloomFilter::trackMemory = true;
		BitVector::trackMemory   = true;
		}
	if (contains(debug,"bvcreation"))
		BitVector::reportCreation = true;

	if (contains(debug,"interval"))
		cerr << "interval is " << startPosition << ".." << endPosition << endl;

	// $$$ consider deterministically shuffling the leaf list, i.e. sort it
	//     .. then shuffle with a seed-based PRNG; this would produce the same
	//     .. tree for the same input regardless of order of the list

	find_leaf_vectors ();

	if (contains(debug,"offsets"))
		{
		for (const auto& bv : leafVectors)
			cerr << "bit vector " << bv->filename << " " << bv->offset << endl;
		}

	// create a binary tree

	cluster_greedily();

	// remove fruitless nodes

	if (cullNodes)
		{
		compute_det_ratio(treeRoot,/*isRoot*/true);
		if (deriveCullingThreshold)
			determine_culling_threshold(treeRoot,/*isRoot*/true);
		cull_nodes(treeRoot,/*isRoot*/true);
		}

	// assign nodes top-down numbers;  nodes will be assigned names from these

	if (renumberNodes)
		top_down_numbering(treeRoot,/*depth*/0,/*isRoot*/true);

	// output the topology

	if (contains(debug,"console"))
		print_topology(cout,treeRoot,0);
	else
		{
	    std::ofstream out(treeFilename);
		print_topology(out,treeRoot,0);
		}

	// clean up

	FileManager::close_file();	// make sure the last bloom filter file we
								// .. opened for read gets closed

	// build the tree (we defer this to the "build" command)

	string commandLine = "howdesbt build \"" + treeFilename + "\"";

	if (inhibitBuild)
		{
		cerr << treeFilename << " has been created"
		     << ", but the internal nodes have not been built." << endl;
		cerr << "You can use this command to build them:" << endl;
		cerr << commandLine << endl;
		}
	else
		deferredCommands.emplace_back(commandLine);

	return EXIT_SUCCESS;
	}

//----------
//
// find_leaf_vectors--
//	Determine the bit vectors that will be the leaves of the tree.
//
// We don't *load* the vectors, but establish a list of BitVector objects that
// point to the subset interval within the corresponding bloom filter file.
//
//----------

void ClusterCommand::find_leaf_vectors()
	{
	// read the filter headers; validate that they all have bit vectors of the
	// correct type, and that they all have the same parameters

	std::ifstream in (listFilename);
	if (not in)
		fatal ("error: failed to open \"" + listFilename + "\"");

	BloomFilter* firstBf = nullptr;
	bool listIsEmpty = true;

	string bfFilename;
	int lineNum = 0;
	while (std::getline (in, bfFilename))
		{
		lineNum++;

		// read the filter's header and verify filter consistency and vector
		// types; note that this does *not* load the bit vector
		// note bene: we keep the first filter open (until we're done) so we
		// can check that all the other vectors are consistent with it
		// $$$ MULTI_VECTOR what if the filter contains more than one bit vector!

		BloomFilter* bf = new BloomFilter (strip_blank_ends(bfFilename));
		bf->preload();
		BitVector* bv = bf->get_bit_vector();
		if (bv->compressor() != bvcomp_uncompressed)
			fatal ("error: bit vectors in \"" + bfFilename + "\" are not uncompressed");

		if (firstBf == nullptr)
			{
			firstBf = bf;
			listIsEmpty = false;
			if (bf->numBits <= startPosition)
				fatal ("error: " + bfFilename + " has only " + std::to_string(bf->numBits) + " bits"
				     + ", so the bit interval "
				     + std::to_string(startPosition) + ".."  + std::to_string(endPosition)
				     + " would be empty");
			if (bf->numBits < endPosition)
				{
				endPosition = bf->numBits;
				cerr << "warning: reducing bit interval to " << startPosition << ".." << endPosition << endl;
				}
			}
		else
			bf->is_consistent_with (firstBf, /*beFatal*/ true);

		// discard the bloom filter (and its bit vector) and create a new "raw"
		// bit vector with the desired bit subset interval

		string bvFilename = bv->filename;
		size_t offset = bv->offset;
		if (bf != firstBf) delete bf;

		size_t startOffset = offset + sdslbitvectorHeaderBytes + startPosition/8;
		string rawFilename = bvFilename
		                   + ":raw"
		                   + ":" + std::to_string(startOffset)
		                   + ":" + std::to_string(endPosition-startPosition);

		BitVector* rawBv = BitVector::bit_vector(rawFilename);
		leafVectors.emplace_back(rawBv);
		}

	if (firstBf != nullptr) delete firstBf;

	in.close();

	// make sure the list wasn't empty

	if (listIsEmpty)
		fatal ("error: \"" + listFilename + "\" contains no bloom filters");
	}

//----------
//
// cluster_greedily--
//	Determine a binary tree structure by greedy clustering.
//
// The clustering process consists of repeatedly (a) choosing the closest pair
// of nodes, and (b) replacing those nodes with a new node that is their union.
//
//----------
//
// Implementation notes:
//	(1)	We use a STL priority queue to keep track of node-to-node distances.
//		Note that by using greater-than as the comparison operator means the
//		smallest distances are on the top of the queue.
//	(2)	We involve the subtree height in the comparison, as a tie breaker. This
//		is to prevent a known degenerate case where a batch of empty nodes (all
//		of which have distance zero to each other) would cluster like a ladder.
//		In a more general case it may keep the overall tree height shorter, but
//		such cases are probably rare.
// 
//----------

struct MergeCandidate
	{
public:
	u64 d;		// distance between u and v
	u32 height;	// height of a subtree containing the u-v merger as its root
	u32 u;		// one node (index into node)
	u32 v;		// other node (index into node)
	};

bool operator>(const MergeCandidate& lhs, const MergeCandidate& rhs)
	{
	if (lhs.d      != rhs.d)      return (lhs.d      > rhs.d);
	if (lhs.height != rhs.height) return (lhs.height > rhs.height);
	if (lhs.u      != rhs.u)      return (lhs.u      > rhs.u);
	return (lhs.v > rhs.v);
	}

void ClusterCommand::cluster_greedily()
	{
	u64 numBits = endPosition - startPosition;
	u64 numBytes = (numBits + 7) / 8;
	u32 numLeaves = leafVectors.size();

	if (numLeaves == 0)
		fatal ("internal error: cluster_greedily() asked to cluster an empty nodelist");

	if (numLeaves == 1)  // (this test is because we will assume the root is not a leaf)
		fatal ("internal error: cluster_greedily() asked to cluster a single node");

	u32 numNodes = 2*numLeaves - 1;  // nodes in tree, including leaves
	BinaryTree* node[numNodes];

	// load the bit arrays for the leaves

	for (u32 u=0 ; u<numLeaves ; u++)
		{
		BitVector* bv = leafVectors[u];
		bv->load();
		node[u] = new BinaryTree(u,bv->bits->data());
		if (node[u] == nullptr)
			fatal ("error: failed to create BinaryTree for node[" + std::to_string(u) + "]");
		if (trackMemory)
			cerr << "@+" << node[u] << " creating BinaryTree for node[" << u << "] (" << bv->filename << ")" << endl;
		node[u]->trackMemory = trackMemory;

		if (contains(debug,"bits"))
			{ cerr << u << ": ";  dump_bits (cerr, node[u]->bits);  cerr << endl; }
		}

	// fill the priority queue with all-vs-all distances among the leaves

	std::priority_queue<MergeCandidate, vector<MergeCandidate>, std::greater<MergeCandidate>> q;

	for (u32 u=0 ; u<numLeaves-1 ; u++)
		{
		for (u32 v=u+1 ; v<numLeaves ; v++)
			{
			u64 d = hamming_distance (node[u]->bits, node[v]->bits, numBits);
			if (contains(debug,"distances"))
				cerr << "node " << u << " vs " << "node " << v << " d=" << d << " h=" << 2 << endl;
			if (contains(debug,"queue"))
				cerr << "pushing (" << d << "," << 2 << "," << u << "," << v << ")" << endl;
			MergeCandidate c = { d,2,u,v };
			q.push (c);
			}
		}

	// for each new node,
	//	- pop the closest active pair (u,v) from the queue
	//	- create a new node w = union of (u,v)
	//	- deactivate u and v by removing their bit arrays
	//	- add the distance to w from each active node

	for (u32 w=numLeaves ; w<numNodes ; w++)
		{
		// pop the closest active pair (u,v) from the queue; nodes that have
		// been deactivate have a null bit array, but still have entries in
		// the queue

		u64 d;
		u32 height, u, v;

		while (true)
			{
			if (q.empty())
				fatal ("internal error: cluster_greedily() queue is empty");

			MergeCandidate cand = q.top();
			q.pop();
			if (contains(debug,"queue"))
				cerr << "popping (" << cand.d << "," << cand.height << "," << cand.u << "," << cand.v << ")"
				     << " q.size()=" << q.size() << endl;
			if (node[cand.u]->bits == nullptr) continue; // u isn't active
			if (node[cand.v]->bits == nullptr) continue; // v isn't active
			d      = cand.d;
			height = cand.height;
			u      = cand.u;
			v      = cand.v;
			break;
			}

		if (contains(debug,"mergings"))
			cerr << "merge " << u << " and " << v << " to make " << w
			     << " (hamming distance " << d << ")" << endl;

		// create a new node w = union of (u,v)

		u64* wBits = (u64*) new char[numBytes];
		if (wBits == nullptr)
			fatal ("error: failed to allocate " + std::to_string(numBytes) + " bytes"
			     + " for node " + std::to_string(w) + "'s bit array");
		if (trackMemory)
			cerr << "@+" << wBits << " allocating bits for node[" << w << "]"
			     << " (merges node[" << u << "] and node[" << v << "])" << endl;

		bitwise_or (node[u]->bits, node[v]->bits, /*dst*/ wBits, numBits);

		node[w] = new BinaryTree(w,wBits,node[u],node[v]);
		if (node[w] == nullptr)
			fatal ("error: failed to create BinaryTree for node[" + std::to_string(w) + "]");
		if (trackMemory)
			cerr << "@+" << node[w] << " creating BinaryTree for node[" << w << "]" << endl;
		node[w]->trackMemory = trackMemory;

		if (contains(debug,"bits"))
			{ cerr << w << ": ";  dump_bits (cerr, wBits);  cerr << endl; }

		// deactivate u and v by removing their bit arrays; if either was a
		// leaf tell the corresonding bit vector it can get rid of its bits;
		//
		// note that if we're going to be culling, we move (or copy) the bit
		// arrays to bCup rather than get rid of them

		if (cullNodes)
			{
			if (u < numLeaves)
				{
				node[u]->bCup = (u64*) new char[numBytes];
				if (node[u]->bCup == nullptr)
					fatal ("error: failed to allocate " + std::to_string(numBytes) + " bytes"
					     + " for node " + std::to_string(u) + "'s bCup array");
				if (trackMemory)
					cerr << "@+" << node[u]->bCup << " allocating bCup for node[" << u << "]" << endl;
				std::memcpy (/*to*/ node[u]->bCup, /*from*/ node[u]->bits, /*how much*/ numBytes);
				}
			else
				{ node[u]->bCup = node[u]->bits; }

			if (v < numLeaves)
				{
				node[v]->bCup = (u64*) new char[numBytes];
				if (node[v]->bCup == nullptr)
					fatal ("error: failed to allocate " + std::to_string(numBytes) + " bytes"
					     + " for node " + std::to_string(v) + "'s bCup array");
				if (trackMemory)
					cerr << "@+" << node[v]->bCup << " allocating bits for node[" << v << "]" << endl;
				std::memcpy (/*to*/ node[v]->bCup, /*from*/ node[v]->bits, /*how much*/ numBytes);
				}
			else
				{ node[v]->bCup = node[v]->bits; }
			}

		if (u < numLeaves)
			leafVectors[u]->discard_bits();
		else if (!cullNodes)
			{
			if (trackMemory)
				cerr << "@-" << node[u]->bits << " discarding bits for node[" << u << "]" << endl;
			delete[] node[u]->bits;
			}

		if (v < numLeaves)
			leafVectors[v]->discard_bits();
		else if (!cullNodes)
			{
			if (trackMemory)
				cerr << "@-" << node[v]->bits << " discarding bits for node[" << v << "]" << endl;
			delete[] node[v]->bits;
			}

		node[u]->bits = nullptr;
		node[v]->bits = nullptr;

		// add the distance to w from each active node

		for (u32 x=0 ; x<w ; x++)
			{
			if (node[x]->bits == nullptr) continue; // x isn't active
			u64 d = hamming_distance (node[x]->bits, wBits, numBits);
			u32 h = 1 + std::max (height,node[x]->height);
			if (contains(debug,"distances"))
				cerr << "node " << x << " vs " << "node " << w << " d=" << d << " h=" << h << endl;
			if (contains(debug,"queue"))
				cerr << "pushing (" << d << "," << h << "," << x << "," << w << ")" << endl;
			MergeCandidate c = { d,h,x,w };
			q.push (c);
			}
		}

	// get rid of the root
	//
	// note that if we're going to be culling, we move (or copy) the bit
	// array to bCup rather than get rid of it

	u32 root = numNodes-1;

	if (cullNodes)
		{  // (note that we assume the root cannot be a leaf)
		node[root]->bCup = node[root]->bits;
		}
	else
		{
		if (trackMemory)
			cerr << "@-" << node[root]->bits << " discarding bits for node[" << root << "]" << endl;
		delete[] node[root]->bits;
		}

	node[root]->bits = nullptr;

	// sanity check -- the only thing left in node list should be the root

	bool failure = false;
	for (u32 x=0 ; x<numNodes ; x++)
		{
		if (node[x]->bits == nullptr) continue; // x isn't active
		cerr << "uh-oh: node " << x << " was never merged" << endl;
		failure = true;
		}

	if (failure)
		fatal ("internal error: cluster_greedily() sanity check failed");

	treeRoot = node[root];
	}

//----------
//
// compute_det_ratio--
//	Collect statistics describing the 'active det ratio'-- the node-by-node
//	fraction of determined-active bits that are determined.
//
//----------
//
// Implementation notes:
//	(1)	The concept of "determined" bits is the same as is used for
//		DeterminedFilter.  But the implementation here shares no code with
//		that, instead making use of the simpler formula
//		  bDet = bCap union complement of bCup
// 
//----------

void ClusterCommand::compute_det_ratio
   (BinaryTree*	node,
	bool		isRoot)
	{
	u64 numBits = endPosition - startPosition;
	u64 numBytes = (numBits + 7) / 8;

	bool isLeaf = (node->children[0] == nullptr);
	if ((node->children[0] == nullptr) != (node->children[1] == nullptr))
		fatal ("internal error: node[" + std::to_string(node->nodeNum) + "] has only one child");

	if (node->bCup == nullptr)
		fatal ("internal error: leaf node[" + std::to_string(node->nodeNum) + "] has no bCup");

	// allocate bit vectors for bCap and bDet

	node->bCap = (u64*) new char[numBytes];
	if (node->bCap == nullptr)
		fatal ("error: failed to allocate " + std::to_string(numBytes) + " bytes"
		     + " for node " + std::to_string(node->nodeNum) + "'s bCap array");
	if (trackMemory)
		cerr << "@+" << node->bCap << " allocating bCap for node[" << node->nodeNum << "]" << endl;

	node->bDet = (u64*) new char[numBytes];
	if (node->bDet == nullptr)
		fatal ("error: failed to allocate " + std::to_string(numBytes) + " bytes"
			 + " for node " + std::to_string(node->nodeNum) + "'s bDet array");
	if (trackMemory)
		cerr << "@+" << node->bDet << " allocating bDet for node[" << node->nodeNum << "]" << endl;

	// if this is a leaf, just copy bCup to bCap, and 'compute' bDet from that
	// $$$ we don't really need bDet, since it will be all ones

	if (isLeaf)
		{
		std::memcpy (/*to*/ node->bCap, /*from*/ node->bCup, /*how much*/ numBytes);
		bitwise_or_not(/*from*/     node->bCap,
					   /*or not*/   node->bCup,
					   /*to*/       node->bDet,
					   /*how much*/ numBits);
		return;
		}

	// otherwise, this is a non-leaf node;  first, process the descendents

	compute_det_ratio(node->children[0]);
	compute_det_ratio(node->children[1]);

	// compute bCap from the children
	//
	// for this node n with children c1 and c2,
	//   bCap(n) = bCap(c0) intersect bCap(c1)

	bitwise_and(/*from*/     node->children[0]->bCap,
				/*and*/      node->children[1]->bCap,
				/*to*/       node->bCap,
				/*how much*/ numBits);

	// compute bDet from bCap and bCup
	//
	// for this node n,
	//   bDet(n) = bCap(n) union (not bCup(n))

	bitwise_or_not(/*from*/     node->bCap,
				   /*or not*/   node->bCup,
				   /*to*/       node->bDet,
				   /*how much*/ numBits);

	// compute det_ratio of the children
	//
	// for each child c,
	//   bDetAct(c) = not bDet(n)   (active bits of bDet at c)
	//   det_ratio = #bDet(c) / #bDetAct(c)
	//             = #(bDet(c) and not bDetAct(c)) / #(not bDet(n))
	//             = #(bDet(c) and not bDet(n))    / (numBits - #bDet(n))

	u64 numDetInf = numBits - bitwise_count(node->bDet,numBits);
	for (int childIx=0 ; childIx<2 ; childIx++)
		{
		BinaryTree* child = node->children[childIx];

		child->numDetOne = bitwise_mask_count(/*in*/         child->bDet,
		                                      /*but not in*/ node->bDet,
		                                      /*how much*/   numBits);
		child->numDetInf = numDetInf;

		if (contains(debug,"detratio"))
			{
			bool childIsLeaf = (child->children[0] == nullptr);
			if (childIsLeaf)
				cerr << "detRatio node[" << child->nodeNum << "]";
			else
				cerr << "detRatio node[" << child->nodeNum << "]"
				     << " (=" << child->children[0]->nodeNum << "+" << child->children[1]->nodeNum << ")";
			cerr << " " << child->numDetOne << "/" << child->numDetInf
				 << " (" << (float(child->numDetOne)/child->numDetInf) << ")"
				 << endl;
			}
		}

	// if this node has no parent, compute its det_ratio
	//
	// det_ratio = #bDet / #bDetAct
	//           = #bDet / numBits

	if (isRoot)
		{
		node->numDetOne = bitwise_count(node->bDet,numBits);
		node->numDetInf = numBits;

		if (contains(debug,"detratio"))
			{
			cerr << "detRatio node[" << node->nodeNum << "]"
			     << " (=" << node->children[0]->nodeNum << "+" << node->children[1]->nodeNum << ")"
			     << " " << node->numDetOne << "/" << node->numDetInf
			     << " (" << (float(node->numDetOne)/node->numDetInf) << ")"
			     << endl;
			}
		}

	// dispose of childrens' bit vectors

	if (not isLeaf)
		{
		for (int childIx=0 ; childIx<2 ; childIx++)
			{
			BinaryTree* child = node->children[childIx];

			if (trackMemory)
				{
				if (child->bCup != nullptr)
					cerr << "@-" << child->bCup << " discarding bCup for node[" << child->nodeNum << "]" << endl;
				if (child->bCap != nullptr)
					cerr << "@-" << child->bCap << " discarding bCap for node[" << child->nodeNum << "]" << endl;
				if (child->bDet != nullptr)
					cerr << "@-" << child->bDet << " discarding bDet for node[" << child->nodeNum << "]" << endl;
				}

			if (child->bCup != nullptr)
				{ delete[] child->bCup;  child->bCup = nullptr; }

			if (child->bCap != nullptr)
				{ delete[] child->bCap;  child->bCap = nullptr; }

			if (child->bDet != nullptr)
				{ delete[] child->bDet;  child->bDet = nullptr; }
			}
		}

	// if this node has no parent, dispose of its bit vectors

	if (isRoot)
		{
		if (trackMemory)
			{
			if (node->bCup != nullptr)
				cerr << "@-" << node->bCup << " discarding bCup for node[" << node->nodeNum << "]" << endl;
			if (node->bCap != nullptr)
				cerr << "@-" << node->bCap << " discarding bCap for node[" << node->nodeNum << "]" << endl;
			if (node->bDet != nullptr)
				cerr << "@-" << node->bDet << " discarding bDet for node[" << node->nodeNum << "]" << endl;
			}

		if (node->bCup != nullptr)
			{ delete[] node->bCup;  node->bCup = nullptr; }

		if (node->bCap != nullptr)
			{ delete[] node->bCap;  node->bCap = nullptr; }

		if (node->bDet != nullptr)
			{ delete[] node->bDet;  node->bDet = nullptr; }
		}

	}

//----------
//
// determine_culling_threshold--
//	Derive a culling threshold from the distribution of active det ratio.
//
// Note that the active det ratio at each is expected to have been computed by
// compute_det_ratio().
// 
//----------

void ClusterCommand::determine_culling_threshold
   (BinaryTree*	node,
	bool		isRoot)
	{
	bool isLeaf = (node->children[0] == nullptr);
	if ((node->children[0] == nullptr) != (node->children[1] == nullptr))
		fatal ("internal error: node[" + std::to_string(node->nodeNum) + "] has only one child");

	// initialize sums at the root

	if (isRoot)
		{
		detRatioSum   = detRatioSumofSquare = 0.0;
		detRatioDenom = 0;
		}

	// add this node's det_ratio to the sums; note that leaves don't contribute
	// to the distribution we're interested in

	if ((not isLeaf) and (node->numDetInf > 0))
		{
		double detRatio = double(node->numDetOne) / node->numDetInf;
		detRatioSum         += detRatio;
		detRatioSumofSquare += detRatio*detRatio;
		detRatioDenom++;

		if (contains(debug,"detratiodistrib"))
			cerr << "detRatio node[" << node->nodeNum << "] " << detRatio << endl;
		}

	// process the descendents

	if (not isLeaf)
		{
		determine_culling_threshold(node->children[0]);
		determine_culling_threshold(node->children[1]);
		}

	// if we're the root, compute the threshold

	if (isRoot)
		{
		if (detRatioDenom == 0)
			fatal ("internal error: can't compute culling threshold, tree has no participating nodes");
		double detRatioMean = detRatioSum/detRatioDenom;
		double detRatioStd  = sqrt(detRatioSumofSquare/detRatioDenom - detRatioMean*detRatioMean);
		cullingThreshold = detRatioMean - cullingThresholdSD*detRatioStd;

		if (contains(debug,"detratiodistrib"))
			cerr << "detRatio mean="  << detRatioMean
			     << " stdev=" << detRatioStd
			     << " cull=" << cullingThreshold
			     << " (across " << detRatioDenom << " nodes)"
			     << endl;

		if (cullingThreshold < 0.0)
			cullingThreshold = 0.0;
		else if (cullingThreshold > 1.0)
			cullingThreshold = 1.0;

		cout << "setting culling threshold to "
		     << std::setprecision(1) << std::fixed << 100*cullingThreshold << "%"
		     << std::setprecision(6) << std::fixed
		     << " (mean=" << detRatioMean << " stdev=" << detRatioStd << ")"
		     << endl;
		}
	}

//----------
//
// cull_nodes--
//	Remove "fruitless" nodes from the clustered binary tree structure.
//
// The culling process consists of removing nodes that have a low percentage of
// bits that will "determine" present/absent for their subtree.  Experiments
// (using real queries) have shown that this ratio correlates well with the
// probability that a node will resolve a query.
//
// Note that the active det ratio at each is expected to have been computed by
// compute_det_ratio(), and cullingThreshold has either been set manually or
// computed by determine_culling_threshold().
// 
//----------
//
// Implementation notes:
//	(1)	Fruitless nodes are left in the tree, but are marked as fruitless so
//		that later operations (such as print_topology) can skip them.
// 
//----------

void ClusterCommand::cull_nodes
   (BinaryTree*	node,
	bool		isRoot)
	{
	bool isLeaf = (node->children[0] == nullptr);
	if ((node->children[0] == nullptr) != (node->children[1] == nullptr))
		fatal ("internal error: node[" + std::to_string(node->nodeNum) + "] has only one child");

	// if this is a leaf, ignore it; leaves are always considered fruitful

	if (isLeaf) return;

	// otherwise, this is a non-leaf node;  first, cull the descendents

	cull_nodes(node->children[0]);
	cull_nodes(node->children[1]);

	// determine whether this node is fruitful
	//   fruitfulness ratio = #bDet(n) / #bDetAct(n)

	if (node->numDetOne < node->numDetInf*cullingThreshold)
		{ // (node->numDetOne/node->numDetInf < cullingThreshold)
		node->fruitful = false;
		if (contains(debug,"cull"))
			cerr << "culling removes node[" << node->nodeNum << "] "
				 << node->numDetOne << "/" << node->numDetInf
				 << " (" << (float(node->numDetOne)/node->numDetInf) << ")"
				 << endl;
		}
	}

//----------
//
// top_down_numbering--
//	Assign node numbers sequentially, starting with the root and moving
//	down level by level, and numbering left to right within each level.
//
//----------

void ClusterCommand::top_down_numbering
   (BinaryTree*	node,
	int			depth,
	bool		isRoot)
	{
	if (isRoot)
		{
		count_depths(node,0);
		int maxDepth = depthToNodeCount.size()-1;
		u32 sum = 0;
		for (int depthIx=0 ; depthIx<=maxDepth ; depthIx++)
			{
			depthToNodeId.emplace_back(sum);
			sum += depthToNodeCount[depthIx];
			}
		}

	bool isLeaf = (node->children[0] == nullptr);
	if (isLeaf) return;

	if (node->fruitful)
		node->nodeId = ++depthToNodeId[depth];

	int nextDepth = depth;
	if (node->fruitful) nextDepth++;

	if (node->children[0] != nullptr) top_down_numbering(node->children[0],nextDepth);
	if (node->children[1] != nullptr) top_down_numbering(node->children[1],nextDepth);
	}

//----------
//
// count_depths--
//	Count the number of nodes at each level of the tree, ignoring leaves and
//	fruitless nodes.
//
// This populates the vector depthToNodeCount, which is assumed to be empty
// upon entry.
//
//----------

void ClusterCommand::count_depths
   (BinaryTree*	node,
	int			depth)
	{
	bool isLeaf = (node->children[0] == nullptr);
	if (isLeaf) return;

	if (node->fruitful)
		{
		while (depthToNodeCount.size() <= (size_t) depth)
			depthToNodeCount.emplace_back(0);
		depthToNodeCount[depth]++;
		}

	int nextDepth = depth;
	if (node->fruitful) nextDepth++;

	if (node->children[0] != nullptr) count_depths(node->children[0],nextDepth);
	if (node->children[1] != nullptr) count_depths(node->children[1],nextDepth);
	}

//----------
//
// print_topology--
//
//----------

void ClusterCommand::print_topology
   (std::ostream&	out,
	BinaryTree*		node,
	int				level)
	{
	u32				numLeaves = leafVectors.size();
	u32				nodeNum = node->nodeNum;
	string			nodeName;

	if (node->fruitful)
		{
		if (nodeNum < numLeaves)
			nodeName = leafVectors[nodeNum]->filename;
		else
			{
			u32 nodeId = node->nodeId;
			if (!renumberNodes) nodeId = 1+nodeNum;

			nodeName = nodeTemplate;
			bool countFromZero = false;
			string field = "{number}";
			std::size_t fieldIx = nodeName.find(field);
			if (fieldIx == string::npos)
				{
				field = "{number:1}";
				fieldIx = nodeName.find(field);
				}
			if (fieldIx == string::npos)
				{
				field = "{number:0}";
				fieldIx = nodeName.find(field);
				countFromZero = (fieldIx != string::npos);
				}

			if (fieldIx == string::npos)
				fatal ("internal error: nodeTemplate=\"" + nodeTemplate + "\""
					 + "does not contain \"{number}\", nor a variant of it");

			if (countFromZero) nodeId--;
			nodeName.replace (fieldIx, field.length(), std::to_string(nodeId));
			}

		if (not contains(debug,"numbers"))
			out << string(level,'*');
		else if (level == 0)
			out << "- (" << nodeNum << ") ";
		else
			out << string(level,'*') << " (" << nodeNum << ") ";
		out << nodeName << endl;
		}

	int nextLevel = level;
	if (node->fruitful) nextLevel++;

	if (node->children[0] != nullptr) print_topology (out, node->children[0], nextLevel);
	if (node->children[1] != nullptr) print_topology (out, node->children[1], nextLevel);
	}

//----------
//
// dump_bits--
//	Write a bit array to a steam, in human-readable form (for debugging).
//
//----------

void ClusterCommand::dump_bits
   (std::ostream&	out,
	const u64*		bits)
	{
	u64		numBits = endPosition - startPosition;
	string	bitsString(numBits,'-');
	u64*	scan = (u64*) bits;

	u64 chunk = *(scan++);
	int bitInChunk = -1;
	for (u64 ix=0 ; ix<numBits ; ix++)
		{
		if (++bitInChunk == 64)
			{ chunk = *(scan++);  bitInChunk = 0; }
		if (((chunk >> bitInChunk) & 1) == 1) bitsString[ix] = '+';
		}

	out << bitsString;
	}
