#ifndef cmd_cluster_H
#define cmd_cluster_H

#include <string>
#include <cstdlib>
#include <cstdint>
#include <iostream>

#include "bit_vector.h"
#include "commands.h"

class BinaryTree
	{
public:
	BinaryTree(const std::uint32_t _nodeNum, std::uint64_t* _bits,
	           BinaryTree* child0=nullptr, BinaryTree* child1=nullptr)
	  :	nodeNum(_nodeNum),
		nodeId(0),
		fruitful(true),
		bits(_bits),
		bCup(nullptr),
		bCap(nullptr),
		bDet(nullptr)
		{
		children[0] = child0;
		children[1] = child1;
		height = 1;
		if (child0 != nullptr) height = 1 + child0->height;
		if (child1 != nullptr) height = std::max(height,1+child1->height);
		}

	virtual ~BinaryTree()
		{
		if (children[0] != nullptr) delete children[0];
		if (children[1] != nullptr) delete children[1];

		if (trackMemory)
			{
			if (bits != nullptr)
				std::cerr << "@-" << bits << " discarding bits for node[" << nodeNum << "]" << std::endl;
			if (bCup != nullptr)
				std::cerr << "@-" << bCup << " discarding bCup for node[" << nodeNum << "]" << std::endl;
			if (bCap != nullptr)
				std::cerr << "@-" << bCap << " discarding bCap for node[" << nodeNum << "]" << std::endl;
			if (bDet != nullptr)
				std::cerr << "@-" << bDet << " discarding bDet for node[" << nodeNum << "]" << std::endl;
			}

		if (bits != nullptr) delete bits;
		if (bCup != nullptr) delete bCup;
		if (bCap != nullptr) delete bCap;
		if (bDet != nullptr) delete bDet;

		if (trackMemory)
			std::cerr << "@-" << this << " discarding BinaryTree node" << std::endl;
		}

	std::uint32_t nodeNum;
	std::uint32_t nodeId;
	bool fruitful;
	std::uint32_t height;
	std::uint64_t* bits;
	std::uint64_t* bCup;		// union of all leaves in the subtree
	std::uint64_t* bCap;		// inersection of all leaves in the subtree
	std::uint64_t* bDet;		// "determined" bits at this node
	BinaryTree* children[2];

	std::uint64_t numDetInf;	// number of active bits in B_det
	std::uint64_t numDetOne;	// number of bits for which B_det==1

	bool trackMemory;
	};


class ClusterCommand: public Command
	{
public:
	static const std::uint64_t defaultEndPosition = 100*1000;
	static constexpr double defaultCullingThresholdSD = 2.0;  // (two standard deviations below the mean)
	//static constexpr double defaultCullingThreshold = 0.2;  // (no longer used)

public:
	ClusterCommand(const std::string& name): Command(name),treeRoot(nullptr) {}
	virtual ~ClusterCommand();
	virtual void short_description (std::ostream& s);
	virtual void usage (std::ostream& s, const std::string& message="");
	virtual void debug_help (std::ostream& s);
	virtual void parse (int _argc, char** _argv);
	virtual int execute (void);
	virtual void find_leaf_vectors (void);
	virtual void cluster_greedily (void);
	virtual void compute_det_ratio (BinaryTree* node,bool isRoot=false);
	virtual void determine_culling_threshold (BinaryTree* node,bool isRoot=false);
	virtual void cull_nodes (BinaryTree* node,bool isRoot=false);
	virtual void top_down_numbering (BinaryTree* node,int depth,bool isRoot=false);
	virtual void count_depths (BinaryTree* node,int depth);
	virtual void print_topology (std::ostream& out, BinaryTree* node, int level);
	virtual void dump_bits (std::ostream& out, const std::uint64_t* bits);

	std::string listFilename;
	std::string treeFilename;
	std::string nodeTemplate;
	std::uint64_t startPosition;	// origin-zero, half-open
	std::uint64_t endPosition;
	bool cullNodes;
	bool deriveCullingThreshold;
	double cullingThresholdSD;
	double cullingThreshold;
	bool renumberNodes;
	bool inhibitBuild;
	bool trackMemory;

	double detRatioSum;
	double detRatioSumofSquare;
	std::uint32_t detRatioDenom;

	std::vector<BitVector*> leafVectors;
	BinaryTree* treeRoot;
	std::vector<std::uint32_t> depthToNodeCount;
	std::vector<std::uint32_t> depthToNodeId;
	};

#endif // cmd_cluster_H
