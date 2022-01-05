// bloom_tree.cc-- classes representing bloom filter trees.
//
// References:
//
//   [1]  Solomon, Brad, and Carl Kingsford. "Fast search of thousands of
//        short-read sequencing experiments." Nature biotechnology 34.3 (2016):
//        300-302.

#include <string>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <iostream>

#include "utilities.h"
#include "bit_utilities.h"
#include "file_manager.h"
#include "bloom_tree.h"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
#define u32 std::uint32_t
#define u64 std::uint64_t

//----------
//
// initialize class variables
//
//----------


//----------
//
// BloomTree--
//
//----------

BloomTree::BloomTree
   (const string& _name,
	const string& _bfFilename)
	  :	isDummy(_bfFilename.empty()),
		manager(nullptr),
		name(_name),
		bfFilename(_bfFilename),
		bf(nullptr),
		isLeaf(true),
		parent(nullptr),
		fpRateKnown(false),
		fpRate(0.0),
		nodesShareFiles(false),
		queryStats(nullptr)
	{
	}

BloomTree::BloomTree
   (BloomTree* root)
	  :	isDummy(root->isDummy),
		manager(nullptr),
		name(root->name),
		bfFilename(root->bfFilename),
		bf(root->bf),
		isLeaf(root->isLeaf),
		parent(nullptr),
		nodesShareFiles(false),
		queryStats(nullptr)
	{
	// nota bene: this doesn't copy the subtree, just the root node; we expect
	//            the caller will detach everything from the root node

	// …… can we do this using copy semantics?
	for (const auto& child : root->children)
		children.emplace_back (child);


	}

BloomTree::~BloomTree()
	{


	if (bf != nullptr) delete bf;
	for (const auto& subtree : children)
		delete subtree;

	if (queryStats != nullptr) delete[] queryStats;
	}

void BloomTree::preload()
	{
	if (bf == nullptr) bf = BloomFilter::bloom_filter(bfFilename);
	bf->preload();
	}

void BloomTree::load()
	{
	if (bf == nullptr)
		{
		bf = BloomFilter::bloom_filter(bfFilename);
		}
	if (manager != nullptr) bf->manager = manager;
	bf->load(/*bypassManager*/false,/*whichNodeName*/name);
	}

void BloomTree::save(bool finished)
	{
	if (bf == nullptr) bf = BloomFilter::bloom_filter(bfFilename);

	for (int bvIx=0 ; bvIx<bf->numBitVectors ; bvIx++)
		{
		BitVector* bv = bf->get_bit_vector(bvIx);
		if (finished)
			{
			bv = bf->simplify_bit_vector(bvIx);
			bv->finished();
			}
		else
			bv->unfinished();
		}

	bf->save();
	}

void BloomTree::unloadable()
	{
	// $$$ eventually we will want a more sophisticated caching mechanism


	if (bf != nullptr)
		{
		if (bf->manager != nullptr)
			bf->discard_bits();
		else
			{ delete bf;  bf = nullptr; }
		}
	}


void BloomTree::add_child
   (BloomTree* offspring)
	{
	children.emplace_back (offspring);
	offspring->parent = this;
	isLeaf = false;
	}

void BloomTree::disown_children()
	{
	children.clear();	// it is assumed the caller has saved references to
						// .. all children prior to asking us to disown them
	}

BloomTree* BloomTree::child
   (size_t	childNum)
	{
	if (childNum >= children.size())
		fatal ("internal error: request for child"
		       " #" + std::to_string(childNum)
		     + " but " + name
		     + " only has " + std::to_string(children.size()) + " children");
	return children[childNum];
	}

BloomFilter* BloomTree::real_filter()
	{
	// skip through dummies to find an instance of a representative bloom
	// filter; for these purposes, we assume all the bloom filters in the tree
	// are similar; note that "real" doesn't mean the filter is loaded or even
	// preloaded

	if (not is_dummy())
		{
		if (bf == nullptr) bf = BloomFilter::bloom_filter(bfFilename);
		if (manager != nullptr) bf->manager = manager;
		return bf;
		}

	if (not isLeaf)
		{
		for (const auto& child : children)
			{
			BloomFilter* realBf = child->real_filter();
			if (realBf != nullptr) return realBf;
			}
		}

	return nullptr;
	}

void BloomTree::pre_order
   (vector<BloomTree*>& order)
	{
	if (not is_dummy()) // (dummies are left out of the resulting list)
		order.emplace_back (this);
	for (const auto& child : children)
		child->pre_order (order);
	}

void BloomTree::post_order
   (vector<BloomTree*>& order)
	{
	for (const auto& child : children)
		child->post_order (order);
	if (not is_dummy()) // (dummies are left out of the resulting list)
		order.emplace_back (this);
	}

void BloomTree::leaves
   (vector<BloomTree*>& order)
	{
	if (isLeaf)
		order.emplace_back (this);
	else
		{
		for (const auto& child : children)
			child->leaves (order);
		}
	}

void BloomTree::print_topology
   (std::ostream&	out,
	int				level,
	int				format) const
	{
	int levelInc = 1;

	if (isDummy)
		levelInc = 0;
	else
		{
		switch (format)
			{
			case topofmt_nodeNames:
				out << string(level,'*') << name << endl;
				break;
			case topofmt_fileNames:
			default:
				out << string(level,'*') << bfFilename << endl;
				break;
			case topofmt_containers:
				out << string(level,'*') << name << "[" << bfFilename << "]" << endl;
				break;
			}
		}

	for (const auto& child : children)
		child->print_topology (out, level+levelInc, format);
	}

//~~~~~~~~~~
// build union tree
//~~~~~~~~~~

void BloomTree::construct_union_nodes (u32 compressor)
	{
	// if we already have a filter, just make sure it is in the pre-loaded
	// state (or beyond)

	if (bf != nullptr)
		{ bf->preload();  return; }

	// if we're compressing, compose a filename for the compressed version of
	// the node;  note that we keep that new name separate from the node's
	// simple name until after we're done with the node;  we write directly to
	// the new name, and by finally installing the new name in the node, the
	// calling program can scan the tree for an accurate topology

	if (compressor != bvcomp_uncompressed)
		{
		string bfKindStr       = "." + BloomFilter::filter_kind_to_string(bfkind_simple);
		string compressionDesc = "." + BitVector::compressor_to_string(compressor);
		if (bfKindStr == ".") bfKindStr = "";
		if (compressionDesc == ".uncompressed") compressionDesc = "";
		futureBfFilename = name + bfKindStr + compressionDesc + ".bf";
		}

	// if this is a leaf, create and load its filter;  if we're NOT compressing
	// we're done;  but if we ARE compressing we write a compressed copy to the
	// new file

	if (isLeaf)
		{

		bf = BloomFilter::bloom_filter(bfFilename);
		bf->load();

		if (compressor != bvcomp_uncompressed)
			{
			if (bf->numBitVectors!=1)
				fatal ("error: " + bfFilename + " contains more than one bit vector");
			BitVector* bvInput = bf->get_bit_vector(0);
			if (bvInput->is_compressed())
				fatal ("error: " + bfFilename + " contains a compressed bit vector");

			BloomFilter newBf(bf,futureBfFilename);
			newBf.new_bits(bvInput,compressor);
			
			newBf.save();
			}

		return;
		}

	// otherwise this is an internal node; first construct its descendants


	for (const auto& child : children)
			child->construct_union_nodes(compressor);

	// if this is a dummy node, we don't need to build it, but we do mark its
	// children as unloadable
	// nota bene: we don't expect a dummy to be a child of some other node

	if (is_dummy())
		{
		for (const auto& child : children)
			{
			child->unloadable();
			if (not child->futureBfFilename.empty())
				{
				child->bfFilename = child->futureBfFilename;
				child->futureBfFilename = "";
				}
			}
		return;
		}

	// create this filter from the union of the child filters, then mark the
	// children as unloadable

	if (bf != nullptr)
		fatal ("internal error: unexpected non-null filter for " + bfFilename);

	bool isFirstChild = true;
	for (const auto& child : children)
		{
		child->load();  // nota bene: child should have already been loaded

		if (child->bf == nullptr)
			fatal ("internal error: failed to load " + child->bfFilename);
		BitVector* childBv = child->bf->get_bit_vector(0);
		if (childBv == nullptr)
			fatal ("internal error: failed to load bit vector for " + child->bfFilename);
		if (childBv->compressor() != bvcomp_uncompressed)
			fatal ("error: " + child->bfFilename + " contains compressed bit vector(s)");


		if (isFirstChild) // incorporate first child's filters
			{
			bf = BloomFilter::bloom_filter(child->bf,bfFilename);
			bf->new_bits(childBv);
			isFirstChild = false;
			}
		else // union with later child's filter
			{
			child->bf->is_consistent_with (bf, /*beFatal*/ true);
			bf->union_with(childBv);
			}

		child->unloadable();
		if (not child->futureBfFilename.empty())
			{
			child->bfFilename = child->futureBfFilename;
			child->futureBfFilename = "";
			}
		}

	if (bf == nullptr)
		fatal ("internal error:"
		       " in construct_union_nodes(\"" + name + "\")"
		     + ", non-leaf node has no children");

	// save the node;  if we're compressing we write a compressed copy to the
	// new file

	if (compressor == bvcomp_uncompressed)
		{
		save(/*finished*/ true);
		}
	else
		{
		BitVector* bvInput = bf->get_bit_vector(0);
		BloomFilter newBf(bf,futureBfFilename);
		newBf.new_bits(bvInput,compressor);
		newBf.save();
		}

	if (parent == nullptr)
		{
		unloadable();
		if (not futureBfFilename.empty())
			{
			bfFilename = futureBfFilename;
			futureBfFilename = "";
			}
		}
	}

//~~~~~~~~~~
// build allsome tree
//~~~~~~~~~~

void BloomTree::construct_allsome_nodes (u32 compressor)
	{
	// if we already have a filter, just make sure it is in the pre-loaded
	// state (or beyond)

	if (bf != nullptr)
		{ bf->preload();  return; }

	string bfKindStr       = "." + BloomFilter::filter_kind_to_string(bfkind_allsome);
	string compressionDesc = "." + BitVector::compressor_to_string(compressor);
	if (compressionDesc == ".uncompressed") compressionDesc = "";
	string newBfFilename = name + bfKindStr + compressionDesc + ".bf";

	// if this is a leaf, create and load its filter
	//   bvs[0] = B'all(x) = B(x)
	//   bvs[1] = B'some(x) = all zeros
	// Note that both of these will be modified when the parent is constructed

	if (isLeaf)
		{

		BloomFilter* bfInput = BloomFilter::bloom_filter(bfFilename);
		bfInput->load();

		if (bfInput->numBitVectors!=1)
			fatal ("error: " + bfFilename + " contains more than one bit vector");
		BitVector* bvInput = bfInput->get_bit_vector(0);
		if (bvInput->is_compressed())
			fatal ("error: " + bfFilename + " contains a compressed bit vector");

		bf = new AllSomeFilter(newBfFilename);
		bf->copy_properties(bfInput);
		bf->steal_bits(bfInput,/*src*/0,/*dst*/0,compressor);
		delete bfInput;

		bf->new_bits(bvcomp_zeros,1);

		// if this leaf has no parent (i.e. it's an orphan), we need to finish
		// it now, the same way we do (later in this function) for any other
		// parentless node

		bool finished = false;
		if ((parent == nullptr) || (parent->is_dummy()))
			finished = true;

		// $$$ we don't necessarily want to save it;  we want to mark it to be
		//     .. saved, but keep it around if we have enough room, since we'll
		//     .. need it to compute its parent, at which time we'll change it

		bfFilename = newBfFilename;
		save(finished);
		unloadable();
		return;
		}

	// otherwise this is an internal node; first construct its descendants


	for (const auto& child : children)
		child->construct_allsome_nodes(compressor);

	// if this is a dummy node, we don't need to build it, but we do mark its
	// children as unloadable
	// nota bene: we don't expect a dummy to be a child of some other node

	if (is_dummy())
		{
		for (const auto& child : children)
			child->unloadable();
		return;
		}

	// create this filter from its child filters

	if (bf != nullptr)
		fatal ("internal error: unexpected non-null filter for " + bfFilename);

	bool isFirstChild = true;
	for (const auto& child : children)
		{
		child->load();

		if (child->bf == nullptr)
			fatal ("internal error: failed to load " + child->bfFilename);
		BitVector* childBvAll  = child->bf->get_bit_vector(0);
		BitVector* childBvSome = child->bf->get_bit_vector(1);
		if ((childBvAll == nullptr) or (childBvSome == nullptr))
			fatal ("internal error: failed to load bit vector(s) for " + child->bfFilename);
		if (childBvAll->is_compressed())
			fatal ("error: " + child->bfFilename + " contains compressed bit vector(s)");
		if ((childBvSome->is_compressed())
		 && (childBvSome->compressor() != bvcomp_zeros)
		 && (childBvSome->compressor() != bvcomp_ones))
			fatal ("error: " + child->bfFilename + " contains compressed bit vector(s)");


		if (isFirstChild) // incorporate first child's filters
			{
			// bvs[0] = Bcap(x) = B'all(child)
			// bvs[1] = Bcup(x) = B'all(child) union B'some(child)
			bf = new AllSomeFilter(child->bf,newBfFilename);
			bf->new_bits(childBvAll,compressor,0);
			bf->new_bits(childBvAll,compressor,1);
			bf->union_with(childBvSome,1);
			isFirstChild = false;
			}
		else // incorporate later child's filters
			{
			// bvs[0] = Bcap(x) = Bcap(x) intersect B'all(child)
			// bvs[1] = Bcup(x) = Bcup(x) union B'all(child) union B'some(child)
			bf->intersect_with(childBvAll,0);
			bf->union_with(childBvAll,1);
			bf->union_with(childBvSome,1);
			}
		}

	if (bf == nullptr)
		fatal ("internal error:"
		       " in construct_allsome_nodes(\"" + name + "\")"
		     + ", non-leaf node has no children");

	// convert this node from Bcap,Bcup to B'all,B'some
	//   bvs[0] = B'all(x)  = Bcap(x), no modification needed
	//   bvs[1] = B'some(x) = Bcup(x) \ Bcap(x)

	BitVector* bvCap = bf->get_bit_vector(0);
	bf->mask_with(bvCap,1);

	// finish the child nodes
	//   bvs[0] = Bsome(c) = B'some(c), no modification needed
	//   bvs[1] = Ball(c)  = B'all(c) \ B'all(x)

	BitVector* bvAll = bf->get_bit_vector(0);
	for (const auto& child : children)
		{
		child->load();

		child->bf->mask_with(bvAll,0);

		child->save(/*finished*/ true);
		child->unloadable();
		}


	// if this node has no parent, we need to finish it now

	bool finished = false;
	if ((parent == nullptr) || (parent->is_dummy()))
		finished = true;

	// $$$ we don't necessarily want to save it;  we want to mark it to be
	//     .. saved, but keep it around if we have enough room, since we'll
	//     .. need it to compute its parent, at which time we'll change it

	bfFilename = newBfFilename;
	save(finished);
	unloadable();
	}

//~~~~~~~~~~
// build determined tree
//~~~~~~~~~~

void BloomTree::construct_determined_nodes (u32 compressor)
	{
	// if we already have a filter, just make sure it is in the pre-loaded
	// state (or beyond)

	if (bf != nullptr)
		{ bf->preload();  return; }

	string bfKindStr       = "." + BloomFilter::filter_kind_to_string(bfkind_determined);
	string compressionDesc = "." + BitVector::compressor_to_string(compressor);
	if (compressionDesc == ".uncompressed") compressionDesc = "";
	string newBfFilename = name + bfKindStr + compressionDesc + ".bf";

	// if this is a leaf, create and load its filter
	//   bvs[0] = Bdet(x) = all ones
	//   bvs[1] = Bhow(x) = B(x)
	// Note that both of these will be modified when the parent is constructed

	if (isLeaf)
		{

		BloomFilter* bfInput = BloomFilter::bloom_filter(bfFilename);
		bfInput->load();

		if (bfInput->numBitVectors!=1)
			fatal ("error: " + bfFilename + " contains more than one bit vector");
		BitVector* bvInput = bfInput->get_bit_vector(0);
		if (bvInput->is_compressed())
			fatal ("error: " + bfFilename + " contains a compressed bit vector");

		bf = new DeterminedFilter(newBfFilename);
		bf->copy_properties(bfInput);
		bf->steal_bits(bfInput,/*src*/0,/*dst*/1,compressor);
		delete bfInput;

		bf->new_bits(compressor,0);
		BitVector* bvDet = bf->get_bit_vector(0);
		bvDet->fill(1);

		// if this leaf has no parent (i.e. it's an orphan), we need to finish
		// it now, the same way we do (later in this function) for any other
		// parentless node

		bool finished = false;
		if ((parent == nullptr) || (parent->is_dummy()))
			{
			BitVector* bvDet = bf->get_bit_vector(0);
			bf->intersect_with(bvDet,1);
			finished = true;
			}

		// $$$ we don't necessarily want to save it;  we want to mark it to be
		//     .. saved, but keep it around if we have enough room, since we'll
		//     .. need it to compute its parent, at which time we'll change it

		bfFilename = newBfFilename;
		save(finished);
		unloadable();
		return;
		}

	// otherwise this is an internal node; first construct its descendants

	for (const auto& child : children)
		child->construct_determined_nodes(compressor);

	// if this is a dummy node, we don't need to build it, but we do mark its
	// children as unloadable
	// nota bene: we don't expect a dummy to be a child of some other node

	if (is_dummy())
		{
		for (const auto& child : children)
			child->unloadable();
		return;
		}

	// create this filter from its child filters
	//   bvs[0] = Bdet(x) = Bcap(x) union complement of Bcup(x)
	//                    = Bhow(x) union z
	//   bvs[1] = Bhow(x) = Bcap(x)
	//                    = intersection, over children c, of Bhow(c)
	//            z       = intersection, over children c, of (Bdet(c) intersect complement of Bhow(c))


	if (bf != nullptr)
		fatal ("internal error: unexpected non-null filter for " + bfFilename);

	bool isFirstChild = true;
	for (const auto& child : children)
		{
		child->load();

		if (child->bf == nullptr)
			fatal ("internal error: failed to load " + child->bfFilename);
		BitVector* childBvDet = child->bf->get_bit_vector(0);
		BitVector* childBvHow = child->bf->get_bit_vector(1);
		if ((childBvHow == nullptr) or (childBvDet == nullptr))
			fatal ("internal error: failed to load bit vector(s) for " + child->bfFilename);
		if ((childBvHow->is_compressed()) || (childBvDet->is_compressed()))
			fatal ("error: " + child->bfFilename + " contains compressed bit vector(s)");


		if (isFirstChild) // incorporate first child's filters
			{
			// bvs[0] = z       = Bdet(c) intersect complement of Bhow(c)
			// bvs[1] = Bhow(x) = Bhow(c)
			bf = new DeterminedFilter(child->bf,newBfFilename);
			bf->new_bits(childBvDet,compressor,0);
			bf->intersect_with_complement(childBvHow,0);
			bf->new_bits(childBvHow,compressor,1);
			isFirstChild = false;
			}
		else // incorporate later child's filters
			{
			// bvs[0] = z       = z intersect Bdet(c) intersect complement of Bhow(c)
			// bvs[1] = Bhow(x) = Bhow(x) intersect Bhow(c)
			bf->intersect_with(childBvDet,0);
			bf->intersect_with_complement(childBvHow,0);
			bf->intersect_with(childBvHow,1);
			}
		}

	if (bf == nullptr)
		fatal ("internal error:"
		       " in construct_determined_nodes(\"" + name + "\")"
		     + ", non-leaf node has no children");

	// convert this node from the temporary vectors computed in the loop
	//   bvs[0] = Bdet(x) = Bhow(x) union z
	//   bvs[1] = Bhow(x), no modification needed

	BitVector* bvHow = bf->get_bit_vector(1);
	bf->union_with(bvHow,0);

	// incorporate bits from this filter, to finish the child nodes
	//   Idet(c) = active bits of Bdet(c)            = complement of Bdet(x)
	//   Ihow(c) = active bits of Bhow(c)            = Bdet(c) intersect Idet(c)
	//   bvs[0]  = Bdet(c) with inactive bits zeroed = Bdet(c) intersect Idet(c)
	//                                               = Bdet(c) intersect complement of Bdet(x)
	//   bvs[1]  = Bhow(c) with inactive bits zeroed = Bhow(c) intersect Ihow(c)
	//                                               = Bhow(c) intersect Bdet(c) intersect Idet(c)
	//                                               = Bhow(c) intersect bvs[0]

	
	BitVector* bvDet = bf->get_bit_vector(0);
	for (const auto& child : children)
		{
		child->load();

		child->bf->intersect_with_complement(bvDet,0);
		BitVector* bvs0 = child->bf->get_bit_vector(0);
		child->bf->intersect_with(bvs0,1);

		child->save(/*finished*/ true);
		child->unloadable();
		}
		

	// if this node has no parent, we need to finish it now
	//   Idet(x) = active bits of Bdet(x)            = all 1s
	//   Ihow(x) = active bits of Bhow(x)            = Bdet(x) intersect Idet(x)
	//                                               = Bdet(x)
	//   bvs[0]  = Bdet(x) with inactive bits zeroed = Bdet(x) intersect Idet(x)
	//                                               = Bdet(x), no modification needed
	//   bvs[1]  = Bhow(x) with inactive bits zeroed = Bhow(x) intersect Ihow(x)
	//                                               = Bhow(x) intersect Bdet(x)

	bool finished = false;
	if ((parent == nullptr) || (parent->is_dummy()))
		{
		BitVector* bvDet = bf->get_bit_vector(0);
		bf->intersect_with(bvDet,1);
		finished = true;
		}

	// $$$ we don't necessarily want to save it;  we want to mark it to be
	//     .. saved, but keep it around if we have enough room, since we'll
	//     .. need it to compute its parent, at which time we'll change it

	bfFilename = newBfFilename;
	save(finished);
	unloadable();
	}

//~~~~~~~~~~
// build determined,brief tree
//~~~~~~~~~~

void BloomTree::construct_determined_brief_nodes (u32 compressor)
	{
	// if we already have a filter, just make sure it is in the pre-loaded
	// state (or beyond)

	if (bf != nullptr)
		{ bf->preload();  return; }

	string bfKindStr       = "." + BloomFilter::filter_kind_to_string(bfkind_determined_brief);
	string compressionDesc = "." + BitVector::compressor_to_string(compressor);
	if (compressionDesc == ".uncompressed") compressionDesc = "";
	string newBfFilename = name + bfKindStr + compressionDesc + ".bf";

	// if this is a leaf, create and load its filter
	//   bvs[0] = Bdet(x) = all ones
	//   bvs[1] = Bhow(x) = B(x)
	// Note that both of these will be modified when the parent is constructed

	if (isLeaf)
		{

		BloomFilter* bfInput = BloomFilter::bloom_filter(bfFilename);
		bfInput->load();

		if (bfInput->numBitVectors!=1)
			fatal ("error: " + bfFilename + " contains more than one bit vector");
		BitVector* bvInput = bfInput->get_bit_vector(0);
		if (bvInput->is_compressed())
			fatal ("error: " + bfFilename + " contains a compressed bit vector");

		bf = new DeterminedBriefFilter(newBfFilename);
		bf->copy_properties(bfInput);
		bf->steal_bits(bfInput,/*src*/0,/*dst*/1,compressor);
		delete bfInput;

		bf->new_bits(compressor,0);
		BitVector* bvDet = bf->get_bit_vector(0);
		bvDet->fill(1);
		bf->get_bit_vector(0)->filterInfo = DeterminedBriefFilter::notSqueezed;
		bf->get_bit_vector(1)->filterInfo = DeterminedBriefFilter::notSqueezed;

		// if this leaf has no parent (i.e. it's an orphan), we need to finish
		// it now, the same way we do (later in this function) for any other
		// parentless node

		bool finished = false;
		if ((parent == nullptr) || (parent->is_dummy()))
			{
			bf->squeeze_by(bvDet,1);

			bf->get_bit_vector(0)->filterInfo = DeterminedBriefFilter::squeezed;
			bf->get_bit_vector(1)->filterInfo = DeterminedBriefFilter::squeezed;
			finished = true;
			}

		// $$$ we don't necessarily want to save it;  we want to mark it to be
		//     .. saved, but keep it around if we have enough room, since we'll
		//     .. need it to compute its parent, at which time we'll change it

		bfFilename = newBfFilename;
		save(finished);
		unloadable();
		return;
		}

	// otherwise this is an internal node; first construct its descendants


	for (const auto& child : children)
		child->construct_determined_brief_nodes(compressor);

	// if this is a dummy node, we don't need to build it, but we do mark its
	// children as unloadable
	// nota bene: we don't expect a dummy to be a child of some other node

	if (is_dummy())
		{
		for (const auto& child : children)
			child->unloadable();
		return;
		}

	// create this filter from its child filters
	//   bvs[0] = Bdet(x) = Bcap(x) union complement of Bcup(x)
	//                    = Bhow(x) union z
	//   bvs[1] = Bhow(x) = Bcap(x)
	//                    = intersection, over children c, of Bhow(c)
	//            z       = intersection, over children c, of (Bdet(c) intersect complement of Bhow(c))


	if (bf != nullptr)
		fatal ("internal error: unexpected non-null filter for " + bfFilename);

	bool isFirstChild = true;
	for (const auto& child : children)
		{
		child->load();

		if (child->bf == nullptr)
			fatal ("internal error: failed to load " + child->bfFilename);
		BitVector* childBvDet = child->bf->get_bit_vector(0);
		BitVector* childBvHow = child->bf->get_bit_vector(1);
		if ((childBvHow == nullptr) or (childBvDet == nullptr))
			fatal ("internal error: failed to load bit vector(s) for " + child->bfFilename);
		if ((childBvHow->is_compressed()) || (childBvDet->is_compressed()))
			fatal ("error: " + child->bfFilename + " contains compressed bit vector(s)");


		if (isFirstChild) // incorporate first child's filters
			{
			// bvs[0] = z       = Bdet(c) intersect complement of Bhow(c)
			// bvs[1] = Bhow(x) = Bhow(c)
			bf = new DeterminedBriefFilter(child->bf,newBfFilename);
			bf->new_bits(childBvDet,compressor,0);
			bf->intersect_with_complement(childBvHow,0);
			bf->new_bits(childBvHow,compressor,1);
			bf->get_bit_vector(0)->filterInfo = DeterminedBriefFilter::notSqueezed;
			bf->get_bit_vector(1)->filterInfo = DeterminedBriefFilter::notSqueezed;
			isFirstChild = false;
			}
		else // incorporate later child's filters
			{
			// bvs[0] = z       = z intersect Bdet(c) intersect complement of Bhow(c)
			// bvs[1] = Bhow(x) = Bhow(x) intersect Bhow(c)
			bf->intersect_with(childBvDet,0);
			bf->intersect_with_complement(childBvHow,0);
			bf->intersect_with(childBvHow,1);
			}
		}

	if (bf == nullptr)
		fatal ("internal error:"
		       " in construct_determined_brief_nodes(\"" + name + "\")"
		     + ", non-leaf node has no children");

	// convert this node from the temporary vectors computed in the loop
	//   bvs[0] = Bdet(x) = Bhow(x) union z
	//   bvs[1] = Bhow(x), no modification needed

	BitVector* bvHow = bf->get_bit_vector(1);
	bf->union_with(bvHow,0);

	// incorporate bits from this filter, to finish the child nodes
	//   Idet(c) = active bits of Bdet(c)             = complement of Bdet(x)
	//   Ihow(c) = active bits of Bhow(c)             = Bdet(c) intersect Idet(c)
	//   bvs[0]  = Bdet(c) with inactive bits removed = Bdet(c) squeeze Idet(c)
	//   bvs[1]  = Bhow(c) with inactive bits removed = Bhow(c) squeeze Ihow(c)

	sdslbitvector* bDetX = bf->get_bit_vector(0)->bits;
	sdslbitvector* iDetC = new sdslbitvector(*bDetX);
	bitwise_complement (/*dst*/ iDetC->data(), bf->numBits);

	for (const auto& child : children)
		{
		child->load();

		sdslbitvector* bHowC = child->bf->get_bit_vector(0)->bits;
		sdslbitvector* iHowC = new sdslbitvector(*bHowC);
		bitwise_and (/*dst*/ iHowC->data(), /*src*/ iDetC->data(), bf->numBits);

		child->bf->squeeze_by(iDetC,0);
		child->bf->squeeze_by(iHowC,1);

		delete iHowC;

		child->bf->get_bit_vector(0)->filterInfo = DeterminedBriefFilter::squeezed;
		child->bf->get_bit_vector(1)->filterInfo = DeterminedBriefFilter::squeezed;

		child->save(/*finished*/ true);
		child->unloadable();
		}

		delete iDetC;
		

	// if this node has no parent, we need to finish it now
	//   Idet(x) = active bits of Bdet(x)             = all 1s
	//   Ihow(x) = active bits of Bhow(x)             = Bdet(x) intersect Idet(x)
	//                                                = Bdet(x)
	//   bvs[0]  = Bdet(x) with inactive bits removed = Bdet(x) squeeze Idet(x)
	//                                                = Bdet(x), no modification needed
	//   bvs[1]  = Bhow(x) with inactive bits removed = Bhow(x) squeeze Ihow(x)
	//                                                = Bhow(x) squeeze Bdet(x)

	bool finished = false;
	if ((parent == nullptr) || (parent->is_dummy()))
		{
		BitVector* bvDet = bf->get_bit_vector(0);
		bf->squeeze_by(bvDet,1);

		bf->get_bit_vector(0)->filterInfo = DeterminedBriefFilter::squeezed;
		bf->get_bit_vector(1)->filterInfo = DeterminedBriefFilter::squeezed;
		finished = true;
		}

	// $$$ we don't necessarily want to save it;  we want to mark it to be
	//     .. saved, but keep it around if we have enough room, since we'll
	//     .. need it to compute its parent, at which time we'll change it

	bfFilename = newBfFilename;
	save(finished);
	unloadable();
	}

//~~~~~~~~~~
// build intersection tree (to assist in debugging)
//~~~~~~~~~~

void BloomTree::construct_intersection_nodes (u32 compressor)
	{
	if (compressor != bvcomp_uncompressed)
		{
		fatal ("internal error: compression isn't implemented for construct_intersection_nodes()");
		// $$$ implementation would parallel construct_union_nodes; but we
		//     just don't need to perform compression on this type of tree
		}

	// if we already have a filter, just make sure it is in the pre-loaded
	// state (or beyond)

	if (bf != nullptr)
		{ bf->preload();  return; }

	string bfKindStr     = "." + BloomFilter::filter_kind_to_string(bfkind_intersection);
	string newBfFilename = name + bfKindStr + ".bf";

	// if this is a leaf, create and pre-load its filter
	// nota bene: we don't usually expect to be called for leaves, but we
	//            handle it anyway

	if (isLeaf)
		{
		bf = BloomFilter::bloom_filter(bfFilename);
		bf->preload();
		return;
		}

	// otherwise this is an internal node; first construct its descendants;
	// note that we ignore leaf children at this point, to reduce our worst
	// case memory footprint


	for (const auto& child : children)
		{
		if (not child->isLeaf)
			child->construct_intersection_nodes(compressor);
		}

	// if this is a dummy node, we don't need to build it, but we do mark its
	// children as unloadable
	// nota bene: we don't expect a dummy to be a child of some other node

	if (is_dummy())
		{
		for (const auto& child : children)
			child->unloadable();
		return;
		}

	// create the filter from the intersection of the child filters, then mark
	// the children as unloadable



	if (bf != nullptr)
		fatal ("internal error: unexpected non-null filter for " + bfFilename);

	bool isFirstChild = true;
	for (const auto& child : children)
		{
		if (child->isLeaf)
			{
			child->bf = BloomFilter::bloom_filter(child->bfFilename);
			child->bf->preload();
			}

		child->load();

		if (child->bf == nullptr)
			fatal ("internal error: failed to load " + child->bfFilename);
		BitVector* childBv = child->bf->get_bit_vector(0);
		if (childBv == nullptr)
			fatal ("internal error: failed to load bit vector for " + child->bfFilename);
		if (childBv->compressor() != bvcomp_uncompressed)
			fatal ("error: " + child->bfFilename + " contains compressed bit vector(s)");

		if (isFirstChild) // copy first child's filter
			{
			bf = BloomFilter::bloom_filter(child->bf,newBfFilename);
			bf->new_bits(childBv);
			isFirstChild = false;
			}
		else // intersection with later child's filter
			{
			child->bf->is_consistent_with (bf, /*beFatal*/ true);
			bf->intersect_with(childBv);
			}

		child->unloadable();
		}

	if (bf == nullptr)
		fatal ("internal error:"
		       " in construct_intersection_nodes(\"" + name + "\")"
		     + ", non-leaf node has no children");

	save();
	}

//~~~~~~~~~~
// query operations
//~~~~~~~~~~

void BloomTree::batch_query
   (vector<Query*>	queries,
	bool			completeSmerCounts)
	{
	// preload a root, and make sure that a leaf-only operation can work with
	// the type of filter we have

	BloomFilter* bf = real_filter();
	if (bf == nullptr)
		fatal ("internal error: batch_query() unable to locate any bloom filter");
	bf->preload();


	// convert the queries to smers/positions

	for (auto& q : queries)
		{
		q->smerize(bf);
		}

	// make a local copy of the query list (consisting of the same instances)
	// while initializing each query's search details; we need a copy because
	// we'll be reordering the list as we move through the tree

	vector<Query*> localQueries;

	for (auto& q : queries)
		{
		u64 numHashes = q->smerHashes.size();
		if (numHashes == 0)
			{
			cerr << "warning: query \"" << q->name << "\" contains no searchable smers" << endl;
			continue; // (queries with no smers are removed from the search)
			}

		q->numPassed     = 0;
		q->numFailed     = 0;
		q->numHashes  	 = numHashes;
		q->numUnresolved = numHashes;
		q->neededToPass  = ceil (q->threshold * numHashes);
		q->neededToFail  = (numHashes - q->neededToPass) + 1;

		localQueries.emplace_back(q);

		}

	// perform the query


	u64 nbActiveQueries = localQueries.size();
	if (nbActiveQueries > 0)
		perform_batch_query(nbActiveQueries, localQueries, completeSmerCounts);
	}

void BloomTree::perform_batch_query
	(u64			nbActiveQueries,
	vector<Query*>	queries,
	bool			completeSmerCounts)
	{
	u64				nbIncomingQueries = nbActiveQueries;
	u64				qIx;

	// skip through dummy nodes

	if (isDummy)
		{

		for (const auto& child : children)
			child->perform_batch_query(nbActiveQueries,queries,completeSmerCounts);
		return;
		}

	// collect some stats

	if (queryStats != nullptr)
		{
		for (qIx=0 ; qIx<nbIncomingQueries ; qIx++)
			{
			Query* q = queries[qIx];
			queryStats[q->batchIx].examined = true;
			}
		}


	// save query state

	for (qIx=0 ; qIx<nbIncomingQueries ; qIx++)
		{
		Query* q = queries[qIx];
		q->numUnresolvedStack.emplace_back(q->numUnresolved);
		q->numPassedStack.emplace_back(q->numPassed);
		q->numFailedStack.emplace_back(q->numFailed);

		}

	// make sure this node's filter is resident

	load();

	// operate on each query in the batch
	//……… ideally, we'd like to perform this for all siblings, then unload the
	//……… .. siblings, before we descend to the siblings' children

	qIx = 0;
	while (qIx < nbActiveQueries)
		{ // note that nbActiveQueries may change during this loop
		Query* q = queries[qIx];
		bool queryPasses = false;
		bool queryFails  = false;


		u64 positionsToTest = q->numUnresolved;
		u64 posIx = 0;
		while (posIx < positionsToTest)
			{
			// each pass through this loop either increases posIx OR decreases
			// positionsToTest
			//
			// Attribution: the technique of swapping resolved positions to the
			// end of the list was inspired by reference [1]

			u64 hashvalue = q->smerHashes[posIx].first;
			size_t hash_position = q->smerHashes[posIx].second;

			bool posIsResolved = true;
			int resolution = lookup(hashvalue);

			if (resolution == BloomFilter::absent)
				{
				if (++q->numFailed >= q->neededToFail)
					{ queryFails = true;  break; }
				}
			else if (resolution == BloomFilter::present)
				{
				q->pos_present_smers.push_back(hash_position);
				q->numPassed++;
				// if we're NOT computing complete smer counts, we can check
				// whether we've observed enough hits to pass this node early
				if ((not completeSmerCounts) and (q->numPassed >= q->neededToPass))
					{ queryPasses = true;  break; }
				}
			else // if (resolution == BloomFilter::unresolved)
				{
				posIsResolved = false;
				}

			// if hashvalue is resolved, swap it with the end of list, and shorten
			// the list; note that we *don't* increase posIx in this case, as
			// the next iteration needs to process the position we just swapped
			// into that slot

			if ((posIsResolved) and (not isLeaf))
				{
				positionsToTest--;
				std::pair<std::uint64_t,std::size_t> tmp_hash_pos = q->smerHashes[posIx];
				q->smerHashes[posIx] = std::move(q->smerHashes[positionsToTest]);
				q->smerHashes[positionsToTest] = std::move(tmp_hash_pos);
				}

			// otherwise, move on to the next hashvalue 

			else
				posIx++;
			}

		q->numUnresolved = positionsToTest;



		// if the query passes, add it to the list of matches for all leaves in
		// this subtree (nb: this 'subtree' may just be a leaf); note that if
		// we're computing complete smer counts, we have to check whether the
		// node passes here because we avoided that test earlier

		if ((completeSmerCounts) and (isLeaf) and (q->numPassed >= q->neededToPass))
			queryPasses = true;
		

		if (queryPasses)
			query_matches_leaves (q);

		// if the query is resolved, swap it with the end of list, and shorten
		// the list; we leave qIx at the same position, as this now points to
		// the query moved from the end of the list
		//
		// otherwise, just move on to the next query

		if (queryPasses or queryFails)
			{
			nbActiveQueries--;
			queries[qIx] = queries[nbActiveQueries];
			queries[nbActiveQueries] = q;
			}
		else
			{
			qIx++;  // move on to the next query
			}

		// collect some stats

		if (queryStats != nullptr)
			{
			querystats* stats = &queryStats[q->batchIx];

			if (queryPasses) stats->passed = true;
			if (queryFails)  stats->failed = true;
			stats->numPassed     = q->numPassed;
			stats->numFailed     = q->numFailed;
			stats->numUnresolved = q->numUnresolved;

			stats->locallyPassed = stats->numPassed;
			stats->locallyFailed = stats->numFailed;
			if (parent != nullptr)
				{
				querystats* parentStats = &parent->queryStats[q->batchIx];
				stats->locallyPassed -= parentStats->numPassed;
				stats->locallyFailed -= parentStats->numFailed;
				}
			}
		}

	// unless we're going to adjust smers/positions, we don't need this node's
	// filter to be resident any more

	bool isPositionAdjustor = bf->is_position_adjustor();
	if (!isPositionAdjustor) unloadable();

	// sanity check: if we're at a leaf, we should have resolved all queries

	if ((isLeaf) and (nbActiveQueries > 0))
		{
		cerr << "internal error: failed to resolve queries at leaf"
			 << " \"" << bfFilename << "\"" << endl;
		cerr << "unresolved queries:";
		for (qIx=0 ; qIx<nbActiveQueries ; qIx++)
			{
			Query* q = queries[qIx];
			if (qIx == 0) cerr << " " << q->name;
			else cerr << ", " << q->name;
			}
		cerr << endl;
		fatal ();
		}

	// adjust smer/position lists as we move down the tree; for most node types
	// this would be a null operation, but for nodes that use rank/select the
	// position values are modified to reflect the removal of inactive bits in
	// the bloom filters
	if (isPositionAdjustor)
		{
		for (qIx=0 ; qIx < nbActiveQueries ; qIx++)
			{
			Query* q = queries[qIx];
			bf->adjust_positions_in_list(q->smerHashes,q->numUnresolved);
			}
		}

	// pass whatever queries remain down to the subtrees

	if (nbActiveQueries > 0)
		{
		for (const auto& child : children)
			child->perform_batch_query(nbActiveQueries,queries,completeSmerCounts);
		}

	// restore smer/position lists as we move up the tree

	if (isPositionAdjustor)
		{
		for (qIx=0 ; qIx<nbActiveQueries ; qIx++)
			{
			Query* q = queries[qIx];
			bf->restore_positions_in_list(q->smerHashes,q->numUnresolved);
			}
		}

	// if we were adjusting smers/positions, we finally don't need this node's
	// filter to be resident any more

	if (isPositionAdjustor) unloadable();

	// restore query state

	for (qIx=0 ; qIx<nbIncomingQueries ; qIx++)
		{
		Query* q = queries[qIx];
		q->numUnresolved = q->numUnresolvedStack.back();
		q->numUnresolvedStack.pop_back();

		assert (q->pos_present_smers.size() == q->numPassed);
		q->numPassed = q->numPassedStack.back();
		q->numPassedStack.pop_back();

		// remove previously added present hash positions
		q->pos_present_smers.resize(q->numPassed);


		q->numFailed = q->numFailedStack.back();
		q->numFailedStack.pop_back();

	

		}

	}

void BloomTree::query_matches_leaves
   (Query* q)
	{
	if (not isLeaf)
		{
		for (const auto& child : children)
			child->query_matches_leaves (q);
		} 
	else
		{
		q->matches.emplace_back (name);
		q->matchesNumPassed.emplace_back (q->numPassed);
		// PIERRE 16 DEC 2021. 
		// Store position of positive smers
		// In q->pos_present_smers_stack (vector<unordered_set<size_t>>)
		std::unordered_set <size_t> local_presentHashes (q->pos_present_smers.begin(), q->pos_present_smers.end());
		q->pos_present_smers_stack.emplace_back(std::move(local_presentHashes));
		
		}
	}


int BloomTree::lookup
   (const u64 hash_value) const
	{
	int resolution = bf->lookup(hash_value);
	if  (resolution != BloomFilter::unresolved)
		return resolution;
	else if (isLeaf)
		return BloomFilter::present;
	else
		return BloomFilter::unresolved;
	}





//----------
//
// read_topology--
//	Read a tree topology from a file and create the corresponding tree
//	object(s). The topology file format is described below.
//
//----------
//
// Arguments:
//	const string&	filename:	The name of the topology file to read from.
//
// Returns:
//	The tree's root node.
//
//----------
//
// Notes:
//	(1) The input format was inspired by reference [1], but is *not* compatible
//		with it.
//	(2)	The input format consists of one line per node, with the nodes listed
//		in pre-order:
//			- a parent is listed before it's children
//			- all nodes in the subtree of the first child are listed before any
//			  nodes of the second child (and second before third, and so on)
//		A node's file name is preceded by a string of asterisks indicating the
//		depth of that node in the tree. For example,
//			root.bf
//			*child1.bf
//			**child3.bf
//			***child5.bf
//			***child6.bf
//			**child4.bf
//			*child2.bf
//	(3)	Nodes can also be listed as a node name followed by a bracketed
//		filename. This facilitates having many filters stored in the same file.
//		In the example here, siblings are stored in one file, so that they can
//		be loaded into memory at the same time.
//			root.bf
//			*child1.bf
//			**child3[child4.bf]
//			***child5.bf
//			***child6[child5.bf]
//			**child4.bf
//			*child2[child1.bf]
//	(4)	The tree needn't be binary; nodes can have more than two children.
//	(5)	The tree may be a forest, in which case a dummy root node is added,
//		having the forest's trees' roots as its children.
//	(6) If the topology filename contains a path, that path is prepended to any
//		node filenames that don't already contain a path.
//	(7)	Upon completion, the tree contains *only* BloomTree nodes, none of the
//		underlying bloom filters are loaded; in fact, the underlying bloom
//		filter files are not read, and need not even exist.
//
//----------

// parse_topology_line

struct ParsedLine
	{
	bool hasProblem;
	std::size_t level;
	string name;
	string bfFilename;
	bool hasBracketedFilename;
	};

static ParsedLine parse_topology_line
   (const string& line,
	const string& basePath)
	{
	ParsedLine p;
	p.hasProblem = true;
	p.level = line.find_first_not_of("*");
	p.bfFilename = strip_blank_ends (line.substr(p.level));
	p.name = "";
	p.hasBracketedFilename = false;

	// if the line is of the form name[filename], split out the relevant parts
	// $$$ probably should use regular expressions for this

	string::size_type lBrackIx = p.bfFilename.find_first_of("[");
	string::size_type rBrackIx = p.bfFilename.find_last_of("]");
	if ((lBrackIx == string::npos) && (rBrackIx == string::npos))
		; // no bracketed expression, do nothing
	else if ((lBrackIx != string::npos) != (rBrackIx != string::npos))
		return p;                               // don't have both brackets
	else // if ((lBrackIx != string::npos) && (rBrackIx != string::npos))
		{
		if ((lBrackIx == 0)                       // empty name
		 or (rBrackIx != p.bfFilename.length()-1) // closing bracket not at end
		 or (rBrackIx == lBrackIx+1))             // empty filename
			return p;
		if (p.bfFilename.find_first_of("[]") != lBrackIx)
			return p;                             // extra bracket before the []
		if (p.bfFilename.find_first_of("[]",lBrackIx+1) != rBrackIx)
			return p;                             // extra bracket between the []
		p.name       = p.bfFilename.substr(0,lBrackIx);
		p.bfFilename = p.bfFilename.substr(lBrackIx+1,rBrackIx-(lBrackIx+1));
		p.hasBracketedFilename = true;
		}

	// if no name was specified, derive one from the filename

	if (p.name == "")
		{
		p.name = strip_file_path(p.bfFilename);
		p.name = BloomFilter::strip_filter_suffix(p.name);
		}

	// if the filename doesn't contain a path, prepend the base path; note that
	// the base path might be empty

	string::size_type slashIx = p.bfFilename.find_last_of("/");
	if (slashIx == string::npos) p.bfFilename = basePath + p.bfFilename;

	p.hasProblem = false;
	return p;
	}

// read_topology--

BloomTree* BloomTree::read_topology
   (const string&	filename)
	{
	std::ifstream in (filename);
	if (not in)
		fatal ("error: failed to open \"" + filename + "\"");

	// extract the base file path (if there is one)

	string basePath("");
	string::size_type slashIx = filename.find_last_of("/");
	if (slashIx != string::npos)
		basePath = filename.substr(0,slashIx+1);

	// create a dummy, filterless, node for the root, whose children will
	// comprise a forest; if the root ends up with a single child, we'll use
	// that child as the root instead

	int numNodes = 0;
	BloomTree* root = new BloomTree();
	bool nodesShareFiles = false;

	// parse the topology file; this first method creates nodes for the entire
	// tree

	vector<BloomTree*> stack;
	stack.push_back (root);

	string line;
	int lineNum = 0;
	while (std::getline (in, line))
		{
		lineNum++;
		line = strip_blank_ends (line);
		if (line.empty()) continue;

		ParsedLine p = parse_topology_line(line, basePath);
		if (p.hasProblem)
			fatal ("error: unable to parse (\"" + filename
					+ "\", line " + std::to_string(lineNum) + ")");

		numNodes++;
		BloomTree* node = new BloomTree(p.name,p.bfFilename);
		if (p.hasBracketedFilename) nodesShareFiles = true;

		if ((numNodes == 0) and (p.level != 0))
			fatal ("error: root must be at level zero (\"" + filename
					+ "\", line " + std::to_string(lineNum) + ")");

		while (stack.size() > p.level+1)
			stack.pop_back ();

		if (p.level+1 != stack.size())
			fatal ("error: tree depth jumps from level " + std::to_string(stack.size()-1)
					+ " to " + std::to_string(p.level+1)
					+ " (\"" + filename + "\", line " + std::to_string(lineNum) + ")");

		BloomTree* parent = stack.back();
		parent->add_child (node);

		stack.push_back (node);
		}
		


	in.close();

	if (root->num_children() == 0)
		{
		fatal ("error: empty tree in \"" + filename + "\"");
		delete root;
		return nullptr;
		}

	// dispose of the dummy root, if it has only one child

	if (root->num_children() == 1)
		{
		BloomTree* oldRoot = root;
		root = root->child(0);
		root->parent = nullptr;
		oldRoot->disown_children();
		delete oldRoot;
		}

	root->nodesShareFiles = nodesShareFiles;
	return root;
	}

