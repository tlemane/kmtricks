#ifndef file_manager_H
#define file_manager_H

#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>

#include "bloom_filter.h"
#include "bloom_tree.h"

//----------
//
// classes in this module--
//
//----------

struct BitVectorInfo
	{
	std::string		name;
	std::uint32_t	compressor;	// compressor identifier (one of bvcomp_xxx)
	size_t			offset;		// offset (into file) of bloom filter's data
	size_t			numBytes;	// number of bytes of data; zero means unknown
	};


class FileManager
	{
public:
	FileManager(BloomTree* root, bool validateConsistency=false);
	virtual ~FileManager();

	virtual void preload_content(const std::string& filename);
	virtual void load_content(const std::string& filename,const std::string& whichNodeName="");


public:
	BloomFilter* modelBf;
	std::unordered_map<std::string,BloomTree*> nameToNode; // hash table
									// .. mapping a node name to the associated
									// .. bloom tree node
	std::unordered_map<std::string,std::vector<std::string>*> filenameToNames;
									// hash table mapping a filename to the
									// .. list of names of nodes to be loaded
									// .. from that file
	std::unordered_map<std::string,bool> alreadyPreloaded; // hash table mapping
									// .. a filename to the to true if the file
									// .. has already been preloaded, false if
									// .. not

public:
	static std::ifstream* open_file  (const std::string& filename,
	                                  std::ios_base::openmode mode=std::ios::in,
	                                  bool positionAtStart=false);
	static void           close_file (std::ifstream* in=nullptr, bool really=false);

public:
	static std::string    openedFilename;
	static std::ifstream* openedFile;
	};

#endif // file_manager_H
