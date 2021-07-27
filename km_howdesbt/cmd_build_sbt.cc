// cmd_build_sbt.cc-- build a sequence bloom tree from a topology file

#include <string>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <vector>
#include <random>

#include "utilities.h"
#include "bloom_tree.h"
#include "file_manager.h"

#include "support.h"
#include "commands.h"
#include "cmd_build_sbt.h"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;


void BuildSBTCommand::short_description
   (std::ostream& s)
	{
	s << commandName << "-- build a sequence bloom tree from a topology file and leaves" << endl;
	}

void BuildSBTCommand::usage
   (std::ostream& s,
	const string& message)
	{
	if (!message.empty())
		{
		s << message << endl;
		s << endl;
		}

	short_description(s);
	s << "usage: " << commandName << " <filename> [options]" << endl;
	//    123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	s << "  <filename>           name of the tree toplogy file" << endl;
	s << "  --outtree=<filename> name of topology file to write tree consisting of the" << endl;
	s << "                       filters built" << endl;
	s << "                       (by default we derive a name for the resulting topology" << endl;
	s << "                       from the input filename; but by default no tree is)" << endl;
	s << "                       written for --simple, as it would be the same as the" << endl;
	s << "                       input tree)" << endl;
	s << "  --simple             create tree nodes as simple bloom filters" << endl;
	s << "                       (this is the default)" << endl;
	s << "  --howde              equivalent to --determined,brief --rrr" << endl;
	s << "  --allsome            create tree nodes as all/some bloom filters" << endl;
	s << "  --determined         create tree nodes as determined/how bloom filters" << endl;
	s << "  --determined,brief   create tree nodes as determined/how, but only store" << endl;
	s << "                       active bits" << endl;
	s << "  --uncompressed       create the nodes as uncompressed bit vector(s)" << endl;
	s << "                       (this is the default)" << endl;
	s << "  --rrr                create the nodes as rrr-compressed bit vector(s)" << endl;
	s << "  --roar               create the nodes as roar-compressed bit vector(s)" << endl;
	}

void BuildSBTCommand::debug_help
   (std::ostream& s)
	{
	s << "--debug= options" << endl;
	s << "  trackmemory" << endl;
	s << "  reportrankselect" << endl;
	s << "  bfsimplify" << endl;
	s << "  btunload" << endl;
	s << "  bfcreation" << endl;
	s << "  bfmanager" << endl;
	s << "  bvcreation" << endl;
	s << "  topology" << endl;
	s << "  load" << endl;
	s << "  traversal" << endl;
	s << "  nochildupdate" << endl;
	}

void BuildSBTCommand::parse
   (int		_argc,
	char**	_argv)
	{
	int		argc;
	char**	argv;

	// defaults

	bfKind     = bfkind_simple;
	compressor = bvcomp_uncompressed;
	BloomTree::inhibitBvSimplify = false;

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

		// --outtree=<filename>

		if (is_prefix_of (arg, "--outtree="))
			{ outTreeFilename = argVal;  continue; }

		// node type

		if ((arg == "--simple")
		 || (arg == "--union")
		 || (arg == "--cup"))
			{ bfKind = bfkind_simple;  continue; }

		if ((arg == "--howde")
		 || (arg == "--HowDe")
		 || (arg == "--howdesbt")
		 || (arg == "--HowDeSBT"))
			{
			bfKind     = bfkind_determined_brief;
			compressor = bvcomp_rrr;
			continue;
			}

		if ((arg == "--allsome")
		 || (arg == "--all/some")
		 || (arg == "--all-some")
		 || (arg == "--all_some"))
			{ bfKind = bfkind_allsome;  continue; }

		if ((arg == "--determined")
		 || (arg == "--determinedhow")
		 || (arg == "--determined/how")
		 || (arg == "--determined-how")
		 || (arg == "--determined_how")
		 || (arg == "--how/de")
		 || (arg == "--how-de")
		 || (arg == "--how_de")
		 || (arg == "--how/det")
		 || (arg == "--how-det")
		 || (arg == "--how_det")
		 || (arg == "--det")
		 || (arg == "--dethow")
		 || (arg == "--det/how")
		 || (arg == "--det-how")
		 || (arg == "--det_how"))
			{ bfKind = bfkind_determined;  continue; }

		if ((arg == "--determined,brief")
		 || (arg == "--determinedhow,brief")
		 || (arg == "--determined/how,brief")
		 || (arg == "--determined-how,brief")
		 || (arg == "--determined_how,brief")
		 || (arg == "--how/de,brief")
		 || (arg == "--how-de,brief")
		 || (arg == "--how_de,brief")
		 || (arg == "--how/det,brief")
		 || (arg == "--how-det,brief")
		 || (arg == "--how_det,brief")
		 || (arg == "--det,brief")
		 || (arg == "--dethow,brief")
		 || (arg == "--det/how,brief")
		 || (arg == "--det-how,brief")
		 || (arg == "--det_how,brief"))
			{ bfKind = bfkind_determined_brief;  continue; }

		// (unadvertised) intersection node type

		if ((arg == "--intersect")
		 || (arg == "--intersection")
		 || (arg == "--cap"))
			{ bfKind = bfkind_intersection;  continue; }

		// compression type

		if (arg == "--uncompressed")
			{ compressor = bvcomp_uncompressed;  continue; }

		if ((arg == "--rrr")
		 || (arg == "--RRR"))
			{ compressor = bvcomp_rrr;  continue; }

		if ((arg == "--roar")
		 || (arg == "--roaring"))
			{ compressor = bvcomp_roar;  continue; }

		// (unadvertised) --tree=<filename>, --topology=<filename>

		if ((is_prefix_of (arg, "--tree="))
		 ||	(is_prefix_of (arg, "--intree="))
		 ||	(is_prefix_of (arg, "--topology=")))
			{
			if (not inTreeFilename.empty())
				chastise ("unrecognized option: \"" + arg + "\""
				          "\ntree topology file was already given as \"" + inTreeFilename + "\"");
			inTreeFilename = argVal;
			continue;
			}

		// (unadvertised) --nobvsimplify

		if (arg == "--nobvsimplify")
			{ BloomTree::inhibitBvSimplify = true;  continue; }

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

		// <filename>

		if (not inTreeFilename.empty())
			chastise ("unrecognized option: \"" + arg + "\""
			          "\ntree topology file was already given as \"" + inTreeFilename + "\"");
		inTreeFilename = arg;
		}

	// sanity checks

	if (inTreeFilename.empty())
		chastise ("you have to provide a tree topology file");

	if ((not outTreeFilename.empty()) and (inTreeFilename.empty()))
		chastise ("cannot use --outtree unless you provide the input tree");

	if (bfKind == bfkind_intersection)
		outTreeFilename = "";
	else if ((bfKind != bfkind_simple) and (outTreeFilename.empty()))
		{
		string bfKindStr = BloomFilter::filter_kind_to_string(bfKind);
		outTreeFilename = strip_file_path(inTreeFilename);
		string::size_type dotIx = outTreeFilename.find_last_of(".");
		if (dotIx == string::npos)
			outTreeFilename = outTreeFilename + "." + bfKindStr + ".sbt";
		else if (is_suffix_of(outTreeFilename,".sbt"))
			outTreeFilename = outTreeFilename.substr(0,dotIx) + "." + bfKindStr + ".sbt";
		else
			outTreeFilename = outTreeFilename + "." + bfKindStr + ".sbt";

		cout << "topology will be written to \"" << outTreeFilename << "\"" << endl;
		}

	return;
	}

int BuildSBTCommand::execute()
	{
	if (contains(debug,"trackmemory"))
		{
		FileManager::trackMemory = true;
		BloomTree::trackMemory   = true;
		BloomFilter::trackMemory = true;
		BitVector::trackMemory   = true;
		}
	if (contains(debug,"reportrankselect"))
		BitVector::reportRankSelect = true;
	if (contains(debug,"bfsimplify"))
		BloomFilter::reportSimplify = true;
	if (contains(debug,"btunload"))
		BloomTree::reportUnload = true;
	if (contains(debug,"bfcreation"))
		BloomFilter::reportCreation = true;
	if (contains(debug,"bfmanager"))
		BloomFilter::reportManager = true;
	if (contains(debug,"bvcreation"))
		BitVector::reportCreation = true;

	BloomTree* root = BloomTree::read_topology(inTreeFilename);

	if (contains(debug,"topology"))
		root->print_topology(cerr,/*level*/0,/*format*/topofmt_nodeNames);

	vector<BloomTree*> order;
	root->post_order(order);
	if (contains(debug,"load"))
		{
		for (const auto& node : order)
			node->reportLoad = true;
		}

	bool hasOnlyChildren = false;
	for (const auto& node : order)
		{
		node->reportSave   = true;
		node->dbgTraversal = (contains(debug,"traversal"));
		// node->dbgInhibitChildUpdate = (contains(debug,"nochildupdate"));

		if (node->num_children() == 1)
			{
			hasOnlyChildren = true;
			cerr << "warning: " << node->child(0)->bfFilename
				<< " is an only child" << endl;
			}
		}
	if (hasOnlyChildren)
		fatal ("error: tree contains at least one only child");

	switch (bfKind)
		{
		case bfkind_simple:
			root->construct_union_nodes (compressor);
			break;
		case bfkind_allsome:
			root->construct_allsome_nodes (compressor);
			break;
		case bfkind_determined:
			root->construct_determined_nodes (compressor);
			break;
		case bfkind_determined_brief:
			root->construct_determined_brief_nodes (compressor);
			break;
		case bfkind_intersection:   // to assist in debugging
			root->construct_intersection_nodes (compressor);
			break;
		default:
			fatal ("error: in BuildSBTCommand::execute():"
			       " bad filter code: \"" + std::to_string(bfKind) + "\"");
		}

	if (not outTreeFilename.empty())
		{
	    std::ofstream out(outTreeFilename);
		root->print_topology(out);
		}

	FileManager::close_file();	// make sure the last bloom filter file we
								// .. opened for read gets closed

	delete root;
	return EXIT_SUCCESS;
	}

