// file_manager.cc-- manage associations between bloom filters and files,
//                   including which filers are resident in memory

#include <string>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <chrono>

#include "utilities.h"
#include "bloom_filter.h"
#include "bloom_tree.h"
#include "file_manager.h"

using std::string;
using std::vector;
using std::pair;
using std::cerr;
using std::endl;
#define u32 std::uint32_t

//----------
//
// initialize class variables
//
//----------


string         FileManager::openedFilename  = "";
std::ifstream* FileManager::openedFile      = nullptr;

//----------
//
// FileManager--
//
//----------

FileManager::FileManager
   (BloomTree*	root,
	bool		validateConsistency)
	  :	modelBf(nullptr)
	{
	// scan the tree, setting up hashes from (a) node name to a node's filter
	// filename and from (b) filter filename to names of nodes within; we also
	// verify that the node names are distinct

	std::unordered_map<string,string> nameToFile;

	vector<BloomTree*> order;
	root->post_order(order);

	for (const auto& node : order)
		node->manager = this;

	for (const auto& node : order)
		{
		if (nameToFile.count(node->name) > 0)
			fatal ("error: tree contains more than one node named \""
			     + node->name 
			     + " (in \"" + node->bfFilename + "\""
			     + " and \"" + nameToFile[node->name] + "\")");

		nameToFile[node->name] = node->bfFilename;
		nameToNode[node->name] = node;

		if (filenameToNames.count(node->bfFilename) == 0)
			{
			filenameToNames[node->bfFilename] = new vector<string>;
			alreadyPreloaded[node->bfFilename] = false;
			}
		filenameToNames[node->bfFilename]->emplace_back(node->name);
		}

	// preload the content headers for every node, file-by-file; this has
	// two side effects -- (1) the bloom filter properties are checked for
	// consistency, and (2) we are installed as every bloom filter's manager

	if (validateConsistency)
		{
		for (auto iter : filenameToNames) 
			{
			string filename = iter.first;
			preload_content (filename);
			}
		}

	}

FileManager::~FileManager()
	{

	if (modelBf != nullptr) delete modelBf;
	for (auto iter : filenameToNames) 
		{
		vector<string>* names = iter.second;

		delete names;
		}
	}

void FileManager::preload_content
   (const string&	filename)
	{
	wall_time_ty startTime;

	if (filenameToNames.count(filename) == 0)
		fatal ("internal error: attempt to preload content from"
		       " unknown file \"" + filename + "\"");

	if (alreadyPreloaded[filename]) return;


	std::ifstream* in = FileManager::open_file(filename,std::ios::binary|std::ios::in,
	                                           /* positionAtStart*/ true);
	if (not *in)
		fatal ("error: FileManager::preload_content()"
			   " failed to open \"" + filename + "\"");
	vector<pair<string,BloomFilter*>> content
		= BloomFilter::identify_content(*in,filename);

	

	vector<string>* nodeNames = filenameToNames[filename];
	if (content.size() != nodeNames->size())
		fatal ("error: \"" + filename + "\""
		     + " contains " + std::to_string(content.size()) + " bloom filters"
		     + ", in conflict with the tree's topology"
		     + " (expected " + std::to_string(nodeNames->size()) + ")");

	for (const auto& templatePair : content)
		{
		string       bfName     = templatePair.first;
		BloomFilter* templateBf = templatePair.second;
		if (not contains (*nodeNames, bfName))
			fatal ("error: \"" + filename + "\""
			     + " contains the bloom filter \"" + bfName + "\""
			     + ", in conflict with the tree's topology");

		BloomTree* node = nameToNode[bfName];
		

		// if the node has already been loaded, leave it be

		if ((node->bf != nullptr) and (node->bf->ready))
			continue;

		// copy the template into the node's filter

		if (node->bf == nullptr)
			{
			node->bf = BloomFilter::bloom_filter(templateBf);
			node->bf->manager = this;
			}
		else // node exists but is not ready
			{
			node->bf->copy_properties(templateBf);
			node->bf->setSizeKnown = templateBf->setSizeKnown;
			node->bf->setSize      = templateBf->setSize;
			}

		node->bf->steal_bits(templateBf);
		delete templateBf;

		// make sure all bloom filters in the tree are consistent

		if (modelBf == nullptr)
			modelBf = BloomFilter::bloom_filter(node->bf);
		else
			node->bf->is_consistent_with(modelBf,/*beFatal*/true);
		}

	alreadyPreloaded[filename] = true;


	FileManager::close_file(in);
	}

void FileManager::load_content
   (const string&	filename,
	const string&	_whichNodeName)
	{
	// ……… when we implement a heap, and dirty bits, we'll need to empty the heap here

	if (filenameToNames.count(filename) == 0)
		fatal ("internal error: attempt to load content from"
		       " unknown file \"" + filename + "\"");

	string whichNodeName = _whichNodeName;
	if (not alreadyPreloaded[filename])
		{
		preload_content(filename);
		whichNodeName = "";  // we will load all nodes in this file
		}

//øøø we only need to load this if it hasn't already been loaded

	vector<string>* nodeNames = filenameToNames[filename];
	for (const auto& nodeName : *nodeNames)
		{
		if ((whichNodeName != "") and (nodeName != whichNodeName))
			continue;
		BloomTree* node = nameToNode[nodeName];
		node->bf->load(/*bypassManager*/ true);
		}

	}

//----------
//
// open_file, close_file--
//	Wrapper for input stream open and close, keeping a file open until some
//	other file is needed.  This (hopefully) saves us the overhead of opening a
//	file, closing it, then opening it again.
//
//----------
//
// CAVEAT:	Any command that reads bit vectors needs to call close_file()
//			before exit, to avoid a memory leak.
//
//----------

std::ifstream* FileManager::open_file
   (const string&			filename,
	std::ios_base::openmode	mode, 
	bool                    positionAtStart)
	{
	if ((openedFile != nullptr) && (filename == openedFilename))
		{
		if (positionAtStart)
			openedFile->seekg(0,openedFile->beg);
		return openedFile;
		}

	if (openedFile != nullptr)
		{
		openedFile->close();
		delete openedFile;
		}

	openedFilename = filename;
	openedFile     = new std::ifstream(filename,mode);
	return openedFile;
	}

void FileManager::close_file
   (std::ifstream*	in,
	bool			really)
	{
	if (openedFile == nullptr) return; // $$$ should this be an error?

	if (in == nullptr)
		{ in = openedFile;  really = true; }
	else if (in != openedFile)
		fatal ("error: FileManager::close_file()"
			   " is asked to close the wrong file");

	if (really)
		{
		openedFile->close();
		delete openedFile;
		openedFilename = "";
		openedFile     = nullptr;
		}
	}

