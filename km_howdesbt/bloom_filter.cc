// bloom_filter.cc-- classes representing bloom filters.
//
// References:
//
//   [1]  Solomon, Brad, and Carl Kingsford. "Improved Search of Large
//        Transcriptomic Sequencing Databases Using Split Sequence Bloom
//        Trees." International Conference on Research in Computational
//        Molecular Biology. Springer, Cham, 2017.
//   [2]  https://en.wikipedia.org/wiki/Bloom_filter#Probability_of_false_positives

#include <string>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <vector>
#include <chrono>

#include "utilities.h"
#include "bit_utilities.h"
#include "hash.h"
#include "file_manager.h"
#include "bloom_filter_file.h"
#include "bloom_filter.h"

using std::string;
using std::vector;
using std::pair;
using std::cerr;
using std::endl;
#define u32 std::uint32_t
#define u64 std::uint64_t

//----------
//
// debugging defines--
//
// In "normal" builds, the triggering defines here are *not* #defined, so that
// no run-time penalty is incurred.
//
//----------

//----------
//
// initialize class variables
//
//----------


//----------
//
// BloomFilter--
//
//----------

BloomFilter::BloomFilter
   (const string& _filename)
	  :	ready(false),
		manager(nullptr),
		filename(_filename),
		smerSize(0),
		hasher1(nullptr),
		hasher2(nullptr),
		numHashes(0),
		hashSeed1(0),
		hashSeed2(0),
		hashModulus(0),
		numBits(0),
		setSizeKnown(false),
		setSize(0),
		numBitVectors(1)
	{
	// nota bene: we clear all maxBitVectors entries (instead of just
	//            numBitVectors), because a subclass using this constructor
	//            may increase numBitVectors
	for (int bvIx=0 ; bvIx<maxBitVectors ; bvIx++) bvs[bvIx] = nullptr;

	}

BloomFilter::BloomFilter
   (const string&	_filename,
	u32				_smerSize,
	u32				_numHashes,
	u64				_hashSeed1,
	u64				_hashSeed2,
	u64				_numBits,
	u64				_hashModulus)
	  :	ready(true),
		manager(nullptr),
		filename(_filename),
		smerSize(_smerSize),
		hasher1(nullptr),
		hasher2(nullptr),
		numHashes(_numHashes),
		hashSeed1(_hashSeed1),
		hashSeed2(_hashSeed2),
		numBits(_numBits),
		setSizeKnown(false),
		setSize(0),
		numBitVectors(1)
	{
	// (see note in first constructor)
	for (int bvIx=0 ; bvIx<maxBitVectors ; bvIx++) bvs[bvIx] = nullptr;

	setup_hashers();
	if (_hashModulus == 0) hashModulus = _numBits;
	                  else hashModulus = _hashModulus;

	}

BloomFilter::BloomFilter
   (const BloomFilter* templateBf,
	const string& newFilename)
	  :	ready(true),
		manager(nullptr),
		smerSize(templateBf->smerSize),
		hasher1(nullptr),
		hasher2(nullptr),
		numHashes(templateBf->numHashes),
		hashSeed1(templateBf->hashSeed1),
		hashSeed2(templateBf->hashSeed2),
		hashModulus(templateBf->hashModulus),
		numBits(templateBf->numBits),
		setSizeKnown(false),
		setSize(0),
		numBitVectors(templateBf->numBitVectors)
	{
	// (see note in first constructor)
	for (int bvIx=0 ; bvIx<maxBitVectors ; bvIx++) bvs[bvIx] = nullptr;

	if (newFilename != "") filename = newFilename;
	                  else filename = templateBf->filename;

	setup_hashers();

	}

BloomFilter::~BloomFilter()
	{

	if (hasher1 != nullptr) delete hasher1;
	if (hasher2 != nullptr) delete hasher2;

	// nota bene: we only consider the first numBitVectors entries; the rest
	//            are never used
	for (int bvIx=0 ; bvIx<numBitVectors ; bvIx++)
		{ if (bvs[bvIx] != nullptr) delete bvs[bvIx]; }
	}

string BloomFilter::identity() const
	{
	return class_identity() + ":\"" + filename + "\"";
	}

void BloomFilter::setup_hashers()
	{
	if ((numHashes > 0) && (hasher1 == nullptr))
		hasher1 = new HashCanonical(smerSize,hashSeed1);
	if ((numHashes > 1) && (hasher2 == nullptr))
		hasher2 = new HashCanonical(smerSize,hashSeed2);
	}

bool BloomFilter::preload(bool bypassManager,bool stopOnMultipleContent)
	{
	// preload usually returns true; a return of false occurs when the
	// file contains more than one BF, and stopOnMultipleContent is true

	if (ready) return true;

	// $$$ should we also write the bvs back to disk?
	for (int bvIx=0 ; bvIx<numBitVectors ; bvIx++)
		{
		if (bvs[bvIx] != nullptr)
			{ delete bvs[bvIx];  bvs[bvIx] = nullptr; }
		}

	if ((manager != nullptr) and (not bypassManager))
		{
		// $$$ this should probably honor stopOnMultipleContent
		manager->preload_content(filename);
		// manager will set ready = true
		}
	else
		{
		wall_time_ty startTime;

		std::ifstream* in = FileManager::open_file(filename,std::ios::binary|std::ios::in,
		                                           /* positionAtStart*/ true);
		if (not *in)
			fatal ("error: " + class_identity() + "::preload()"
				   " failed to open \"" + filename + "\"");


		vector<pair<string,BloomFilter*>> content
		    = BloomFilter::identify_content(*in,filename);
		if ((stopOnMultipleContent) and (content.size() != 1))
			{
			FileManager::close_file(in);
			return false;
			}
		if (content.size() != 1)
			fatal ("(internal?) error: in " + identity() + ".preload()"
			     + " file contains multiple bloom filters"
			     + " but we aren't using a file manager");

		BloomFilter* templateBf = content[0].second;
		u32 bfKind = kind();
		u32 templateBfKind = templateBf->kind();
		if (templateBfKind != bfKind)
			fatal ("(internal?) error: in " + identity() + ".preload()"
			     + " file contains incompatible"
			     + "\n.. bloom filter, expected kind=" + filter_kind_to_string(bfKind,false)
			     + " but file has kind=" + filter_kind_to_string(templateBfKind,false));

		copy_properties(templateBf);
		setSizeKnown = templateBf->setSizeKnown;
		setSize      = templateBf->setSize;
		steal_bits(templateBf);
		delete templateBf;
		FileManager::close_file(in);
		}

	if ((numHashes > 0) && (hasher1 == nullptr))
		setup_hashers();

	return true;
	}

//……… we should have a reload() too
//……… we shouldn't be able to load() a filter that didn't come from a file
void BloomFilter::load
   (bool			bypassManager,
	const string&	whichNodeName)
	{
//…… enable this test, non-null and resident and dirty
//	if (bv != nullptr)
//		fatal ("internal error for " + identity()
//		     + "; attempt to load() onto non-null bit vector");
//……

	if ((manager != nullptr) and (not bypassManager))
		{
		manager->load_content(filename,whichNodeName);
		}
	else
		{
		// $$$ assert that whichNodeName == ""
		if (not ready) preload ();

		for (int bvIx=0 ; bvIx<numBitVectors ; bvIx++)
			{
			BitVector* bv = bvs[bvIx];
			bv->load();
			}
		}
	}

void BloomFilter::save()
	{
	for (int bvIx=0 ; bvIx<numBitVectors ; bvIx++)
		{
		if (bvs[bvIx] == nullptr)
			{
			if (bvIx == 0)
				fatal ("internal error for " + identity()
				     + "; attempt to save null bloom filter");
			else
				fatal ("internal error for " + identity()
				     + "; attempt to save partially null bloom filter");
			}
		}

// $$$ this needs to make sure the file isn't currently opened for read!!!!

	// allocate the header, with enough room for a bfvectorinfo record for each
	// component
	//
	// note that we are assuming that the header size for the current file
	// format version is at least as large as that for any earlier versions

	u64 headerBytesNeeded = bffileheader_size(numBitVectors);
	headerBytesNeeded = round_up_16(headerBytesNeeded);

	u32 headerSize = (u32) headerBytesNeeded;
	if (headerSize != headerBytesNeeded)
		fatal ("error: header record for \"" + filename + "\""
		       " would be too large (" + std::to_string(headerSize) + " bytes)");

	bffileheader* header = (bffileheader*) new char[headerSize];
	if (header == nullptr)
		fatal ("error:"
		       " failed to allocate " + std::to_string(headerSize) + " bytes"
		     + " for header record for \"" + filename + "\"");


	// write a fake header to the file; after we write the rest of the file
	// we'll rewind and write the real header; we do this because we don't know
	// the component offsets and sizes until after we've written them


	memset (header, 0, headerSize);
	header->magic      = bffileheaderMagicUn;
	header->headerSize = (u32) sizeof(bffileprefix);

	std::ofstream out (filename, std::ios::binary | std::ios::trunc | std::ios::out);
	out.write ((char*)header, headerSize);
	size_t bytesWritten = headerSize; // (based on assumption of success)
	if (not out)
		fatal ("error: " + class_identity() + "::save(" + identity() + ")"
		     + " failed to open \"" + filename + "\"");

	// start the real header

	header->magic        = bffileheaderMagic;
	header->headerSize   = headerSize;
	header->version      = bffileheaderVersion;
	header->bfKind       = kind();
	header->padding1     = 0;
	header->smerSize     = smerSize;
	header->numHashes    = numHashes;
	header->hashSeed1    = hashSeed1;
	header->hashSeed2    = hashSeed2;
	header->hashModulus  = hashModulus;
	header->numBits      = numBits;
	header->numVectors   = numBitVectors;
	header->setSizeKnown = setSizeKnown;
	header->setSize      = (setSizeKnown)? setSize : 0;

	// write the component(s)

	for (int bvIx=0 ; bvIx<numBitVectors ; bvIx++)
		{
		BitVector* bv = bvs[bvIx];

		header->info[bvIx].compressor = bv->compressor();
		header->info[bvIx].name       = 0;
		header->info[bvIx].offset     = bytesWritten;
		if ((header->info[bvIx].compressor == bvcomp_rrr)
		 || (header->info[bvIx].compressor == bvcomp_unc_rrr))
			{
			header->info[bvIx].compressor |= (RRR_BLOCK_SIZE  << 8);
			header->info[bvIx].compressor |= (RRR_RANK_PERIOD << 16);
			}

		size_t numBytes = bv->serialized_out (out, filename, header->info[bvIx].offset);
		bytesWritten += numBytes;

		header->info[bvIx].numBytes   = numBytes;
		header->info[bvIx].filterInfo = bv->filterInfo;
		}

	// rewind and overwrite header

	out.seekp(std::ios::beg);
	out.write ((char*)header, headerSize);
	out.close();

	// clean up


	if (header != nullptr) delete[] header;

	// now we're in the equivalent of the "ready" state; actually we're beyond
	// that state, in the same state as the result of load()

	ready = true;
	}

void BloomFilter::copy_properties
   (const BloomFilter*	templateBf)
	{
	smerSize    = templateBf->smerSize;
	numHashes   = templateBf->numHashes;
	hashSeed1   = templateBf->hashSeed1;
	hashSeed2   = templateBf->hashSeed2;
	hashModulus = templateBf->hashModulus;
	numBits     = templateBf->numBits;
	}

void BloomFilter::steal_bits
   (BloomFilter*	templateBf)
	{
	if (numBitVectors != templateBf->numBitVectors)
		fatal ("internal error for " + identity()
			 + "; source filter has "
			 + std::to_string(templateBf->numBitVectors) + " bitvectors"
			 + " (this filter has " + std::to_string(numBitVectors) + ")");

	discard_bits();

	for (int bvIx=0 ; bvIx<numBitVectors ; bvIx++)
		{
		bvs[bvIx] = templateBf->bvs[bvIx];
		templateBf->bvs[bvIx] = nullptr;
		}

	ready = true;
	}

void BloomFilter::steal_bits
   (BloomFilter*	templateBf,
	int				whichBv)
	{
	steal_bits(templateBf,/*src*/whichBv,/*dst*/whichBv);
	}

void BloomFilter::steal_bits
   (BloomFilter*	templateBf,
	int				whichSrcBv,
	int				whichDstBv,
	u32				compressor)
	{
	if ((whichDstBv < 0) || (whichDstBv >= numBitVectors))
		fatal ("internal error for " + identity()
			 + "; request to set bitvector " + std::to_string(whichDstBv));
	if ((whichSrcBv < 0) || (whichSrcBv >= templateBf->numBitVectors))
		fatal ("internal error for " + identity()
			 + "; request to get source filter's bitvector " + std::to_string(whichSrcBv));

	discard_bits(whichDstBv);

	BitVector* srcBv = templateBf->bvs[whichSrcBv];
	templateBf->bvs[whichSrcBv] = nullptr;

	if (compressor == srcBv->compressor())
		bvs[whichDstBv] = srcBv;
	else
		{
		bvs[whichDstBv] = BitVector::bit_vector(compressor,srcBv);
		delete srcBv;
		}

	ready = true;
	}

bool BloomFilter::is_consistent_with
   (const BloomFilter*	bf,
	bool				beFatal) const
	{
	if (bf->smerSize != smerSize)
		{
		if (not beFatal) return false;
		fatal ("error: inconsistent smer size " + std::to_string(bf->smerSize)
			 + " in \"" + bf->filename + "\""
			 + " (expected " + std::to_string(smerSize)
			 + " like in \"" + filename + "\")"
			 + "\n(all bloom filters are required to have the same smer size)");
		}

	if (bf->numHashes != numHashes)
		{
		if (not beFatal) return false;
		fatal ("error: inconsistent number of hashes " + std::to_string(bf->numHashes)
			 + " in \"" + bf->filename + "\""
			 + " (expected " + std::to_string(numHashes)
			 + " like in \"" + filename + "\")"
			 + "\n(all bloom filters are required to have the same number of hashes)");
		}

	if (bf->hashSeed1 != hashSeed1)
		{
		if (not beFatal) return false;
		fatal ("error: inconsistent hash seed " + std::to_string(bf->hashSeed1)
			 + " in \"" + bf->filename + "\""
			 + " (expected " + std::to_string(hashSeed1)
			 + " like in \"" + filename + "\")"
			 + "\n(all bloom filters are required to have the same hash seeds)");
		}

	if (bf->hashSeed2 != hashSeed2)
		{
		if (not beFatal) return false;
		fatal ("error: inconsistent hash seed 2 " + std::to_string(bf->hashSeed2)
			 + " in \"" + bf->filename + "\""
			 + " (expected " + std::to_string(hashSeed2)
			 + " like in \"" + filename + "\")"
			 + "\n(all bloom filters are required to have the same hash seeds)");
		}

	if (bf->hashModulus != hashModulus)
		{
		if (not beFatal) return false;
		fatal ("error: inconsistent hash modulus " + std::to_string(bf->hashModulus)
			 + " in \"" + bf->filename + "\""
			 + " (expected " + std::to_string(hashModulus)
			 + " like in \"" + filename + "\")"
			 + "\n(all bloom filters are required to have the same hash modulus -- the same"
			 + "\nnumber of bits)");
		}

	if (bf->numBits != numBits)
		{
		if (not beFatal) return false;
		fatal ("error: inconsistent number of bits " + std::to_string(bf->numBits)
			 + " in \"" + bf->filename + "\""
			 + " (expected " + std::to_string(numBits)
			 + " like in \"" + filename + "\")"
			 + "\n(all bloom filters are required to have the same number of bits)");
		}

	if (bf->kind() != kind())
		{
		if (not beFatal) return false;
		fatal ("error: inconsistent bloom filter kind " + std::to_string(bf->kind())
			 + " in \"" + bf->filename + "\""
			 + " (expected " + std::to_string(kind())
			 + " like in \"" + filename + "\")"
			 + "\n(all bloom filters are required to be of the same kind)");
		}

	return true;
	}

void BloomFilter::discard_bits()
	{
	for (int bvIx=0 ; bvIx<numBitVectors ; bvIx++)
		{
		if (bvs[bvIx] != nullptr)
			{ delete bvs[bvIx];  bvs[bvIx] = nullptr; }
		}
	}

void BloomFilter::discard_bits
   (int whichBv)
	{
	if ((whichBv < 0) || (whichBv >= numBitVectors))
		fatal ("internal error for " + identity()
		     + "; request to discard bitvector " + std::to_string(whichBv));

	if (bvs[whichBv] != nullptr)
		{ delete bvs[whichBv];  bvs[whichBv] = nullptr; }
	}

void BloomFilter::new_bits
   (u32		compressor,
	int		whichBv)
	{
	if ((whichBv < -1) || (whichBv >= numBitVectors))
		fatal ("internal error for " + identity()
		     + "; request to replace bitvector " + std::to_string(whichBv));

	if (whichBv >= 0)
		{
		if (bvs[whichBv] != nullptr) delete bvs[whichBv];
		bvs[whichBv] = BitVector::bit_vector(compressor,numBits);
		}
	else
		{
		// nothing specified, so do them all
		for (int bvIx=0 ; bvIx<numBitVectors ; bvIx++)
			{
			if (bvs[bvIx] != nullptr) delete bvs[bvIx];
			bvs[bvIx] = BitVector::bit_vector(compressor,numBits);
			}
		}
	}

void BloomFilter::new_bits
   (BitVector*	srcBv,
	u32			compressor,
	int			whichBv)
	{
	if ((whichBv < 0) || (whichBv >= numBitVectors))
		fatal ("internal error for " + identity()
		     + "; request to set bitvector " + std::to_string(whichBv));

	if (bvs[whichBv] != nullptr) delete bvs[whichBv];

	if (srcBv->bits == nullptr)
		{
		u32 srcCompressor = srcBv->compressor();
		if ((srcCompressor != bvcomp_zeros) && (srcCompressor != bvcomp_ones))
			fatal ("internal error for " + identity()
			     + "; attempt to copy bits from null or compressed bitvector " + srcBv->identity());
		}

	bvs[whichBv] = BitVector::bit_vector(compressor,srcBv);
	}

void BloomFilter::new_bits
   (const string& filename)  // note that compressor etc. may be encoded in filename
	{
	for (int bvIx=0 ; bvIx<numBitVectors ; bvIx++)
		{
		if (bvs[bvIx] != nullptr) delete bvs[bvIx];
		bvs[bvIx] = BitVector::bit_vector (filename);
		}
	}

BitVector* BloomFilter::get_bit_vector
   (int whichBv) const
	{
	if ((whichBv < 0) || (whichBv >= numBitVectors))
		fatal ("internal error for " + identity()
		     + "; request to get bitvector " + std::to_string(whichBv));

	return bvs[whichBv];
	}

BitVector* BloomFilter::surrender_bit_vector
   (int	 whichBv)
	{
	if ((whichBv < 0) || (whichBv >= numBitVectors))
		fatal ("internal error for " + identity()
		     + "; request to get bitvector " + std::to_string(whichBv));

	BitVector* bv = bvs[whichBv];
	bvs[whichBv] = nullptr;

	return bv;
	}

BitVector* BloomFilter::simplify_bit_vector
   (int	 whichBv)
	{
	// if possible, replace the bit vector with a simpler version, such as
	// all-zeros or all-ones
	if ((whichBv < 0) || (whichBv >= numBitVectors))
		fatal ("internal error for " + identity()
		     + "; request to simplify bitvector " + std::to_string(whichBv));

	BitVector* bv = bvs[whichBv];
	u32 bvCompressor = bv->compressor();

	if ((bvCompressor == bvcomp_zeros) || (bvCompressor == bvcomp_ones))
		return bv;  // bv is already a simple type

	if (bv->is_all_zeros())
		{
		bvs[whichBv] = new ZerosBitVector(bv->size());
		delete bv;
		return bvs[whichBv];
		}

	if (bv->is_all_ones())
		{
		bvs[whichBv] = new OnesBitVector(bv->size());
		delete bv;
		return bvs[whichBv];
		}

	return bv;
	}

void BloomFilter::complement
   (int			whichDstBv)
	{
	if ((whichDstBv < -1) || (whichDstBv >= numBitVectors))
		fatal ("internal error for " + identity()
		     + "; request to complement bitvector " + std::to_string(whichDstBv));

	if (whichDstBv >= 0)
		bvs[whichDstBv]->complement();
	else // if (whichDstBv == -1)
		{
		for (int bvIx=0 ; bvIx<numBitVectors ; bvIx++)
			bvs[bvIx]->complement();
		}
	}

void BloomFilter::union_with
   (BitVector*	srcBv,
	int			whichDstBv)
	{
	if ((whichDstBv < 0) || (whichDstBv >= numBitVectors))
		fatal ("internal error for " + identity()
		     + "; request to union into bitvector " + std::to_string(whichDstBv));

	switch (srcBv->compressor())
		{
		case bvcomp_zeros:
			break;
		case bvcomp_ones:
			bvs[whichDstBv]->fill(1);
			break;
		default:
			bvs[whichDstBv]->union_with(srcBv->bits);
			break;
		}
	}

void BloomFilter::union_with_complement
   (BitVector*	srcBv,
	int			whichDstBv)
	{
	if ((whichDstBv < 0) || (whichDstBv >= numBitVectors))
		fatal ("internal error for " + identity()
		     + "; request to union into bitvector " + std::to_string(whichDstBv));

	switch (srcBv->compressor())
		{
		case bvcomp_zeros:
			bvs[whichDstBv]->fill(1);
			break;
		case bvcomp_ones:
			break;
		default:
			bvs[whichDstBv]->union_with_complement(srcBv->bits);
			break;
		}
	}

void BloomFilter::intersect_with
   (BitVector*	srcBv,
	int			whichDstBv)
	{
	if ((whichDstBv < 0) || (whichDstBv >= numBitVectors))
		fatal ("internal error for " + identity()
		     + "; request to intersection into bitvector " + std::to_string(whichDstBv));

	switch (srcBv->compressor())
		{
		case bvcomp_zeros:
			bvs[whichDstBv]->fill(0);
			break;
		case bvcomp_ones:
			break;
		default:
			bvs[whichDstBv]->intersect_with(srcBv->bits);
			break;
		}

	}

void BloomFilter::mask_with
   (BitVector*	srcBv,
	int			whichDstBv)
	{
	if ((whichDstBv < 0) || (whichDstBv >= numBitVectors))
		fatal ("internal error for " + identity()
		     + "; request to mask bitvector " + std::to_string(whichDstBv));

	switch (srcBv->compressor())
		{
		case bvcomp_zeros:
			break;
		case bvcomp_ones:
			bvs[whichDstBv]->fill(0);
			break;
		default:
			bvs[whichDstBv]->mask_with(srcBv->bits);
			break;
		}
	}

void BloomFilter::xor_with
   (BitVector*	srcBv,
	int			whichDstBv)
	{
	if ((whichDstBv < 0) || (whichDstBv >= numBitVectors))
		fatal ("internal error for " + identity()
		     + "; request to xor into bitvector " + std::to_string(whichDstBv));

	switch (srcBv->compressor())
		{
		case bvcomp_zeros:
			break;
		case bvcomp_ones:
			bvs[whichDstBv]->complement();
			break;
		default:
			bvs[whichDstBv]->xor_with(srcBv->bits);
			break;
		}
	}

void BloomFilter::squeeze_by
   (BitVector*	srcBv,
	int			whichDstBv)
	{
	if ((whichDstBv < 0) || (whichDstBv >= numBitVectors))
		fatal ("internal error for " + identity()
		     + "; request to squeeze bitvector " + std::to_string(whichDstBv));

	// …… require that whichDstBv is uncompressed?

	u32 compressor = srcBv->compressor();
	switch (compressor)
		{
		default:
			bvs[whichDstBv]->squeeze_by(srcBv->bits);
			break;
		case bvcomp_zeros:
		case bvcomp_ones:
			int fillValue = (compressor == bvcomp_zeros)? 0 : 1;
			u64 resultNumBits = bitwise_count(srcBv->bits->data(),numBits);
			sdslbitvector* resultBits = new sdslbitvector(resultNumBits,fillValue);

			bvs[whichDstBv]->replace_bits(resultBits);
			break;
		}
	}

void BloomFilter::squeeze_by
   (const sdslbitvector* srcBits,
	int			whichDstBv)
	{
	if ((whichDstBv < 0) || (whichDstBv >= numBitVectors))
		fatal ("internal error for " + identity()
		     + "; request to squeeze bitvector " + std::to_string(whichDstBv));

	bvs[whichDstBv]->squeeze_by(srcBits);
	}

// mer_to_hash_value--
//	Report the position of a smer in the filter; we return BloomFilter::npos if
//	the smer's position is not within the filter (this is *not* the same as the
//	smer being present in the set represent by the filter); the smer can be a
//	string or 2-bit encoded data

u64 BloomFilter::mer_to_hash_value
   (const string& mer) const
	{
	// nota bene: we don't enforce numHashes == 1

	u64 hashvalue = hasher1->hash(mer) % hashModulus;
	if (hashvalue < numBits) return hashvalue;
	              else return npos;  // position is *not* in the filter
	}

u64 BloomFilter::mer_to_hash_value
   (const u64* merData) const
	{
	// nota bene: we don't enforce numHashes == 1

	u64 hashvalue = hasher1->hash(merData) % hashModulus;
	if (hashvalue < numBits) return hashvalue;
	              else return npos;  // position is *not* in the filter
	}

// add--
//	Add a smer to the filter; the smer can be a string or 2-bit encoded data

void BloomFilter::add
   (const string& mer)
	{
	// nota bene: we don't enforce mer.length == smerSize
	// nota bene: we don't canonicalize the string; revchash handles that

	BitVector* bv = bvs[0];

	u64 h1  = hasher1->hash(mer);
	u64 hashvalue = h1 % hashModulus;
	if (hashvalue < numBits)
		{
		(*bv).write_bit(hashvalue);
		}

	if (numHashes > 1)
		{
		u64 hashValues[numHashes];
		Hash::fill_hash_values(hashValues,numHashes,h1,hasher2->hash(mer));
		for (u32 h=1 ; h<numHashes ; h++)
			{
			hashvalue = hashValues[h] % hashModulus;
			if (hashvalue < numBits)
				{
				(*bv).write_bit(hashvalue); // $$$ MULTI_VECTOR write each bit to a different vector
				}
			}
		}
	}

void BloomFilter::add
   (const u64* merData)
	{
	BitVector* bv = bvs[0];

	u64 h1  = hasher1->hash(merData);
	u64 hashvalue = h1 % hashModulus;
	if (hashvalue < numBits)
		{
		(*bv).write_bit(hashvalue);
		}

	if (numHashes > 1)
		{
		u64 hashValues[numHashes];
		Hash::fill_hash_values(hashValues,numHashes,h1,hasher2->hash(merData));
		for (u32 h=1 ; h<numHashes ; h++)
			{
			hashvalue = hashValues[h] % hashModulus;
			if (hashvalue < numBits)
				{
				(*bv).write_bit(hashvalue); // $$$ MULTI_VECTOR write each bit to a different vector
				}
			}
		}
	}

// contains--
//	returns true if the bloom filter contains the given smer (or false
//	positive), false otherwise; the smer can be a string or 2-bit encoded data

bool BloomFilter::contains
   (const string& mer) const
	{
	// nota bene: we don't enforce mer.length == smerSize
	// nota bene: we don't canonicalize the string; revchash handles that

	BitVector* bv = bvs[0];

	u64 h1  = hasher1->hash(mer);
	u64 hashvalue = h1 % hashModulus;
	if (hashvalue < numBits)
		{
		if ((*bv)[hashvalue] == 0) return false;
		}

	if (numHashes > 1)
		{
		u64 hashValues[numHashes];
		Hash::fill_hash_values(hashValues,numHashes,h1,hasher2->hash(mer));
		for (u32 h=1 ; h<numHashes ; h++)
			{
			hashvalue = hashValues[h] % hashModulus;
			if (hashvalue < numBits)
				{
				if ((*bv)[hashvalue] == 0) return false; // $$$ MULTI_VECTOR read each bit from a different vector
				}
			}
		}

	return true;
	}

bool BloomFilter::contains
   (const u64* merData) const
	{
	BitVector* bv = bvs[0];

	u64 h1  = hasher1->hash(merData);
	u64 hashvalue = h1 % hashModulus;
	if (hashvalue < numBits)
		{
		if ((*bv)[hashvalue] == 0) return false;
		}

	if (numHashes > 1)
		{
		u64 hashValues[numHashes];
		Hash::fill_hash_values(hashValues,numHashes,h1,hasher2->hash(merData));
		for (u32 h=1 ; h<numHashes ; h++)
			{
			hashvalue = hashValues[h] % hashModulus;
			if (hashvalue < numBits)
				{
				if ((*bv)[hashvalue] == 0) return false; // $$$ MULTI_VECTOR read each bit from a different vector
				}
			}
		}

	return true;
	}

int BloomFilter::lookup
   (const u64 hashvalue) const
	{
	BitVector* bv = bvs[0];

	// we assume, without checking, that 0 <= hashvalue < numBits
	if ((*bv)[hashvalue] == 0) return absent;
	else return unresolved;
	}

//----------
//
// AllSomeFilter--
//
//----------

AllSomeFilter::AllSomeFilter
   (const string& _filename)
	  :	BloomFilter(_filename)
	{
	numBitVectors = 2;
	}

AllSomeFilter::AllSomeFilter
   (const string&	_filename,
	u32				_smerSize,
	u32				_numHashes,
	u64				_hashSeed1,
	u64				_hashSeed2,
	u64				_numBits,
	u64				_hashModulus)
	  :	BloomFilter(_filename, _smerSize,
	                _numHashes, _hashSeed1, _hashSeed2,
	                _numBits, _hashModulus)
	{
	numBitVectors = 2;
	}

AllSomeFilter::AllSomeFilter
   (const BloomFilter* templateBf,
	const string& newFilename)
	  :	BloomFilter(templateBf, newFilename)
	{
	numBitVectors = 2;
	}

AllSomeFilter::~AllSomeFilter()
	{
	}

// add--

void AllSomeFilter::add
   (const string& mer)
	{
	fatal ("internal error: attempt to add a mer to " + class_identity());
	}

void AllSomeFilter::add
   (const u64* merData)
	{
	fatal ("internal error: attempt to add a mer to " + class_identity());
	}

// contains--

bool AllSomeFilter::contains
   (const string& mer) const
	{
	fatal ("internal error: \"is mer contained\" request in " + class_identity());
	return false;
	}

bool AllSomeFilter::contains
   (const u64* merData) const
	{
	fatal ("internal error: \"is mer contained\" request in " + class_identity());
	return false;
	}

int AllSomeFilter::lookup
   (const u64 hashvalue) const
	{
	BitVector* bvAll  = bvs[0];
	BitVector* bvSome = bvs[1];

	// we assume, without checking, that 0 <= hashvalue < numBits
	if      ((*bvAll) [hashvalue] == 1) return present;
	else if ((*bvSome)[hashvalue] == 0) return absent;
	else                          return unresolved;
	}

//----------
//
// DeterminedFilter--
//
//----------

DeterminedFilter::DeterminedFilter
   (const string& _filename)
	  :	BloomFilter(_filename)
	{
	numBitVectors = 2;
	}

DeterminedFilter::DeterminedFilter
   (const string&	_filename,
	u32				_smerSize,
	u32				_numHashes,
	u64				_hashSeed1,
	u64				_hashSeed2,
	u64				_numBits,
	u64				_hashModulus)
	  :	BloomFilter(_filename, _smerSize,
	                _numHashes, _hashSeed1, _hashSeed2,
	                _numBits, _hashModulus)
	{
	numBitVectors = 2;
	}

DeterminedFilter::DeterminedFilter
   (const BloomFilter* templateBf,
	const string& newFilename)
	  :	BloomFilter(templateBf, newFilename)
	{
	numBitVectors = 2;
	}

DeterminedFilter::~DeterminedFilter()
	{
	}

int DeterminedFilter::lookup
   (const u64 hashvalue) const
	{
	BitVector* bvDet = bvs[0];
	BitVector* bvHow = bvs[1];

	// we assume, without checking, that 0 <= hashvalue < numBits
	if      ((*bvDet)[hashvalue] == 0) return unresolved;
	else if ((*bvHow)[hashvalue] == 1) return present;
	else                         return absent;
	}

//----------
//
// DeterminedBriefFilter--
//
// Attribution: the use of rank/select and removal of inactive bits was
// inspired by [1], but the application of it to a determined/brief split is
// original with this program.
//
//----------

DeterminedBriefFilter::DeterminedBriefFilter
   (const string& _filename)
	  :	DeterminedFilter(_filename)
	{
	numBitVectors = 2;
	}

DeterminedBriefFilter::DeterminedBriefFilter
   (const string&	_filename,
	u32				_smerSize,
	u32				_numHashes,
	u64				_hashSeed1,
	u64				_hashSeed2,
	u64				_numBits,
	u64				_hashModulus)
	  :	DeterminedFilter(_filename, _smerSize,
	                _numHashes, _hashSeed1, _hashSeed2,
	                _numBits, _hashModulus)
	{
	numBitVectors = 2;
	}

DeterminedBriefFilter::DeterminedBriefFilter
   (const BloomFilter* templateBf,
	const string& newFilename)
	  :	DeterminedFilter(templateBf, newFilename)
	{
	numBitVectors = 2;
	}

DeterminedBriefFilter::~DeterminedBriefFilter()
	{
	}


int DeterminedBriefFilter::lookup
   (const u64 hashvalue) const
	{
	// we assume, without checking, that 0 <= hashvalue < numBits

	BitVector* bvDet = bvs[0];

	if ((*bvDet)[hashvalue] == 0) return unresolved;

	BitVector* bvHow = bvs[1];
	u64 howPos = bvDet->rank1(hashvalue);

	if ((*bvHow)[howPos] == 1) return present;
	else                       return absent;
	}

void DeterminedBriefFilter::adjust_positions_in_list
   (std::vector<std::pair<std::uint64_t,std::size_t>> &smerHashes,
	u64 numUnresolved)
	{
	BitVector* bvDet = bvs[0];

	for (u64 posIx=0 ; posIx<numUnresolved ; posIx++)
		{
		uint64_t hashvalue  = smerHashes[posIx].first;
		uint64_t rank = bvDet->rank1(hashvalue);
		smerHashes[posIx].first = hashvalue - rank;  // $$$ isn't this just rank0(hashvalue)?
		}
	}

void DeterminedBriefFilter::restore_positions_in_list
   (std::vector<std::pair<std::uint64_t,std::size_t>> &smerHashes,
	u64 numUnresolved)
	{
	BitVector* bvDet = bvs[0];

	for (u64 posIx=0 ; posIx<numUnresolved ; posIx++)
		{
		u64 hashvalue = smerHashes[posIx].first;
		smerHashes[posIx].first = bvDet->select0(hashvalue);
		}
	}

//----------
//
// strip_filter_suffix--
//	Remove any of the standard bloom filter suffixes from a file name.
//
//----------
//
// Arguments:
//	const string&	filename:	The name of the bloom filter file.
//	const int		complete:	How much to remove
//								  1 => remove, e.g. ".bf"
//								  2 => remove, e.g. ".rrr.bf"
//								  3 => remove, e.g. ".det.bf" or ".det.rrr.bf"
//
// Returns:
//	A filename, with the suffix removed.
//
//----------

string BloomFilter::strip_filter_suffix
   (const string&	filename,
	const int		complete)
	{
	string name = filename;

	if (is_suffix_of (name, ".bf"))
		name = strip_suffix(name,".bf");

	if (is_suffix_of (name, ".unity"))
		name = strip_suffix(name,".unity");

	if (complete >= 2)
		{
		string compSuffix = "." + BitVector::compressor_to_string(bvcomp_rrr);
		if (is_suffix_of(name,compSuffix))
			name = strip_suffix(name,compSuffix);
		else
			{
			compSuffix = "." + BitVector::compressor_to_string(bvcomp_roar);
			if (is_suffix_of(name,compSuffix))
				name = strip_suffix(name,compSuffix);
			}
		}

	if (complete >= 3)
		{
		if (is_suffix_of(name,".allsome"))			// bfkind_allsome
			name = strip_suffix(name,".allsome");
		else if (is_suffix_of(name,".det"))			// bfkind_determined
			name = strip_suffix(name,".det");
		else if (is_suffix_of(name,".detbrief"))	// bfkind_determined_brief
			name = strip_suffix(name,".detbrief");
		}

	return name;
	}

//----------
//
// default_filter_name--
//	Derive a name for a bloom filter, from its file name.
//
//----------
//
// Arguments:
//	const string&	filename:			The name of the bloom filter file.
//	const int		componentNumber:	The index of the bloom filter in the
//										file's components.  The value -1
//										indicates that this index is not
//										relevant.
//
// Returns:
//	A name for the filter.
//
//----------

string BloomFilter::default_filter_name
   (const string&	filename,
	const int		componentNumber)
	{
	string name = strip_file_path(filename);
	name = BloomFilter::strip_filter_suffix(name);
	if (componentNumber >= 0)
		name += "." + std::to_string(componentNumber);

	return name;
	}

//----------
//
// filter_kind_to_string--
//	Convert a bloom filter kind to a string, often used for constructing
//	filenames.
//
//----------
//
// Arguments:
//	u32		bfKind:			The type of bloom filter; one of bfkind_xxx.
//	bool	shortString:	true  => produce a short string, suitable for use
//							         .. as a file extension
//							false => produce a longer string, suitable for
//							         .. error text
//
// Returns:
//	A short string representing the bloom filter type.
//
//----------

string BloomFilter::filter_kind_to_string
   (u32		bfKind,
	bool	shortString)
	{
	switch (bfKind)
		{
		case bfkind_simple:           return (shortString)? ""         : "simple";
		case bfkind_allsome:          return (shortString)? "allsome"  : "allsome";
		case bfkind_determined:       return (shortString)? "det"      : "determined";
		case bfkind_determined_brief: return (shortString)? "detbrief" : "determined,brief";
		case bfkind_intersection:     return (shortString)? "cap"      : "intersection";
		// nota bene: when new file extensions are added here, they should also
		// .. be added in BloomFilter::strip_filter_suffix()
		default:
			fatal ("error: in filter_kind_to_string():"
			       " bad filter code: \"" + std::to_string(bfKind) + "\"");
		}

	return "";  // should never get here
	}

//----------
//
// vectors_per_filter--
//	Figure out how many bit vectors a particular type of filter will have.
//
//----------
//
// Arguments:
//	const u32 bfKind:	The type of bloom filter; one of bfkind_xxx.
//
// Returns:
//	The number of bit vectors that a filter of that type will have.  This will
//	be the same as the numBitVectors instance variable for the corresponding
//	filter class.
//
//----------

int BloomFilter::vectors_per_filter
   (const u32 bfKind)
	{
	switch (bfKind)
		{
		case bfkind_simple:
		case bfkind_intersection:
			return 1;
			break;
		case bfkind_allsome:
		case bfkind_determined:
		case bfkind_determined_brief:
			return 2;
			break;
		default:
			fatal ("error: in vectors_per_filter():"
			       " bad filter code: \"" + std::to_string(bfKind) + "\"");
		}

	return 0; // never reaches here
	}

//----------
//
// bloom_filter--
//	Create a specified kind of BloomFilter object.
//
//----------
//
// Arguments (variant 1):
//	const std::string& filename:	A filename for the new bloom filter.  The
//									.. type of the bloom filter is derived from
//									.. the suffix.
//
// Arguments (variant 2):
//	const u32 bfKind:				The type of bloom filter to create; one of
//									.. bfkind_xxx.
//	filename..hashModulus:			The same as for the BloomFilter constructors.
//
// Arguments (variant 3):
//	const BloomFilter* templateBf:	A bloom filter to mimic; only in the sense
//									.. of its type and properties.
//	const std::string& newFilename:	A different filename for the new bloom
//									.. filter.
//
// Returns:
//	A BloomFilter (or subclass of BloomFilter) object.
//
//----------

//=== variant 1 ===

BloomFilter* BloomFilter::bloom_filter
   (const string&	filename)
	{
	string reducedName = filename;
	
	string compSuffix = "." + BitVector::compressor_to_string(bvcomp_rrr) + ".bf";
	if (is_suffix_of(reducedName,compSuffix))
		reducedName = strip_suffix(reducedName,compSuffix) + ".bf";
	else
		{
		compSuffix = "." + BitVector::compressor_to_string(bvcomp_roar) + ".bf";
		if (is_suffix_of(reducedName,compSuffix))
			reducedName = strip_suffix(reducedName,compSuffix) + ".bf";
		}

	if (is_suffix_of(reducedName,".detbrief.bf"))		// bfkind_determined_brief
		return new DeterminedBriefFilter (filename);
	if (is_suffix_of(reducedName,".det.bf"))			// bfkind_determined
		return new DeterminedFilter (filename);
	if (is_suffix_of(reducedName,".allsome.bf"))		// bfkind_allsome
		return new AllSomeFilter (filename);
	if (is_suffix_of(reducedName,".bf"))				// bfkind_simple
		return new BloomFilter (filename);

	fatal ("error: BloomFilter::bloom_filter(\"" + filename + "\""
	     + " is not implemented (file extension not recognized)");

	return nullptr; // never reaches here
	}

//=== variant 2 ===

BloomFilter* BloomFilter::bloom_filter
   (const u32		bfKind,
	const string&	filename,
	u32				smerSize,
	u32				numHashes,
	u64				hashSeed1,
	u64				hashSeed2,
	u64				numBits,
	u64				hashModulus)
	{
	switch (bfKind)
		{
		case bfkind_simple:
		case bfkind_intersection:  // internally treated the same as bfkind_simple
			return new BloomFilter      (filename,smerSize,
			                             numHashes,hashSeed1,hashSeed2,
			                             numBits,hashModulus);
		case bfkind_allsome:
			return new AllSomeFilter    (filename,smerSize,
			                             numHashes,hashSeed1,hashSeed2,
			                             numBits,hashModulus);
		case bfkind_determined:
			return new DeterminedFilter (filename,smerSize,
			                             numHashes,hashSeed1,hashSeed2,
			                             numBits,hashModulus);
		case bfkind_determined_brief:
			return new DeterminedBriefFilter
			                            (filename,smerSize,
			                             numHashes,hashSeed1,hashSeed2,
			                             numBits,hashModulus);
		default:
			fatal ("error: BloomFilter::bloom_filter(\"" + std::to_string(bfKind) + "\""
			     + " is not implemented");
		}

	return nullptr; // never reaches here
	}

//=== variant 3 ===

BloomFilter* BloomFilter::bloom_filter
   (const BloomFilter* templateBf,
	const std::string& newFilename)
	{
	switch (templateBf->kind())
		{
		case bfkind_simple:
		case bfkind_intersection:  // internally treated the same as bfkind_simple
			return new BloomFilter(templateBf,newFilename);
		case bfkind_allsome:
			return new AllSomeFilter(templateBf,newFilename);
		case bfkind_determined:
			return new DeterminedFilter(templateBf,newFilename);
		case bfkind_determined_brief:
			return new DeterminedBriefFilter(templateBf,newFilename);
		default:
			fatal ("error: BloomFilter::bloom_filter()"
			       " doesn't understand filter type "
			     + std::to_string(templateBf->kind()));
		}

	return nullptr; // never reaches here
	}

//----------
//
// identify_content--
//	Read a header from a bloom filter file, and determine the content of the
//	file.  "Content" consists of one-or-more named bloom filters.
//
//----------
//
// Arguments:
//	std::ifstream&	in:			The stream to read from.  This should have been
//								opened with std::ios::binary | std::ios::in.
//	const string&	filename:	The name of the bloom filter file to read from.
//								This is only used for error messages, and can
//								be blank.
//
// Returns:
//	The contents; this is a vector of (name,bloom filter) pairs.  The name is
//	the bloom filter's name from the file header, if the header has one.  If
//	not, it is derived from the file name.
//
//----------
//
// Notes:
//	(1)	Each bloom filter created also has a bit vector  (or bit vectors)
//		created for it.  The bit vector has the proper information about where
//		to find its bits (file, offset, number of bytes, compression type), but
//		does not have its bits loaded.
//
//----------

vector<pair<string,BloomFilter*>> BloomFilter::identify_content
   (std::ifstream&	in,
	const string&	filename)
	{
	wall_time_ty startTime;
	double elapsedTime = 0.0;
	vector<pair<string,BloomFilter*>> content;

	//---- read and validate the header prefix ----

	bffileprefix prefix;

	in.read ((char*) &prefix, sizeof(prefix));

	if (!in.good())
		fatal ("error: BloomFilter::identify_content(" + filename + ")"
		       " problem reading header from \"" + filename + "\"");
	size_t prevFilePos    = 0;
	size_t currentFilePos = in.gcount();
	if (currentFilePos != sizeof(prefix))
		fatal ("error: BloomFilter::identify_content(" + filename + ")"
		       " read(\"" + filename + "\"," + std::to_string(sizeof(prefix)) + ")"
		     + " produced " + std::to_string(currentFilePos) + " bytes");
	

	if (prefix.magic == bffileheaderMagicUn)
		fatal ("error: BloomFilter::identify_content(" + filename + ")"
		       " looks like an incomplete bloom filter file"
		       " (it seems the file was not completely written)");

	if (prefix.magic != bffileheaderMagic)
		fatal ("error: BloomFilter::identify_content(" + filename + ")"
		       " doesn't look like a bloom filter file"
		       " (incorrect magic number)");

	if ((prefix.version != bffileheaderVersion)
	 && (prefix.version != bffileheaderVersion1))
		fatal ("error: BloomFilter::identify_content(" + filename + ")"
		       " bloom filter file version " + std::to_string(prefix.version)
		     + " is not supported by this program");

	if (prefix.headerSize <= sizeof(prefix))
		fatal ("error: BloomFilter::identify_content(" + filename + ")"
		       " header impossibly small (" + std::to_string(prefix.headerSize) + " bytes)");

	if (prefix.headerSize > max_bffileheader_size)
		fatal ("error: BloomFilter::identify_content(" + filename + ")"
		       " headers larger than " + std::to_string(max_bffileheader_size) + " bytes are not supported,"
		     + " this file's header claims to be " + std::to_string(prefix.headerSize) + " bytes");

	//---- read the rest of the header, and validate ----

	bffileheader* header = (bffileheader*) new char[prefix.headerSize];
	if (header == nullptr)
		fatal ("error: BloomFilter::identify_content(" + filename + ")"
			   " failed to allocate " + std::to_string(prefix.headerSize) + " bytes"
			 + " for file header");
	std::memcpy (/*to*/ header, /*from*/ &prefix, /*how much*/ sizeof(prefix));


	size_t remainingBytes = prefix.headerSize - sizeof(prefix);
	in.read (((char*) header) + sizeof(prefix), remainingBytes);
	prevFilePos = currentFilePos;  currentFilePos += in.gcount();
	if (currentFilePos != prefix.headerSize)
		fatal ("error: BloomFilter::identify_content(" + filename + ")"
		       " read(\"" + filename + "\"," + std::to_string(remainingBytes) + ")"
		     + " produced " + std::to_string(currentFilePos-prevFilePos) + " bytes");


	if ((header->bfKind != bfkind_simple)
	 && (header->bfKind != bfkind_allsome)
	 && (header->bfKind != bfkind_determined)
	 && (header->bfKind != bfkind_determined_brief)
	 && (header->bfKind != bfkind_intersection))
		fatal ("error: BloomFilter::identify_content(" + filename + ")"
		       " bad filter type: " + std::to_string(header->bfKind));

	// make sure header size and number of vectors make sense

	size_t minHeaderSize = bffileheader_size(header->numVectors);
	if (header->headerSize < minHeaderSize)
		fatal ("error: BloomFilter::identify_content(" + filename + ")"
		       " expected " + std::to_string(minHeaderSize) + " byte header (or larger)"
		     + " but header says it is " + std::to_string(header->headerSize) + " bytes");

	if (header->numVectors < 1)
		fatal ("error: BloomFilter::identify_content(" + filename + ")"
		       " bad number of vectors: " + std::to_string(header->numVectors));

	int vectorsPerFilter = BloomFilter::vectors_per_filter(header->bfKind);
	int numFilters       = header->numVectors / vectorsPerFilter;
	if (header->numVectors % vectorsPerFilter != 0)
		fatal ("error: BloomFilter::identify_content(" + filename + ")"
		       " number of vectors (" + std::to_string(header->numVectors) + ")"
		     + " is not a multiple of the number of vectors per filter"
		     + "(" + std::to_string(vectorsPerFilter) + ")");

	if (header->padding1 != 0)
		fatal ("error: BloomFilter::identify_content(" + filename + ")"
		       " non-zero padding field 1: " + std::to_string(header->padding1));

	if (prefix.version == bffileheaderVersion1)
		{
		u32* padding2 = (u32*) &header->setSizeKnown;
		u32* padding3 = (u32*) &header->setSize;
		u32* padding4 = (u32*) ((((char*)padding3)+sizeof(u32)));

		if (*padding2 != 0)
			fatal ("error: BloomFilter::identify_content(" + filename + ")"
			       " non-zero padding field 2: " + std::to_string(*padding2));

		if (*padding3 != 0)
			fatal ("error: BloomFilter::identify_content(" + filename + ")"
			       " non-zero padding field 3: " + std::to_string(*padding3));

		if (*padding4 != 0)
			fatal ("error: BloomFilter::identify_content(" + filename + ")"
			       " non-zero padding field 4: " + std::to_string(*padding4));
		}

	if (header->numHashes < 1)
		fatal ("error: BloomFilter::identify_content(" + filename + ")"
		       " bad number of hash functions: " + std::to_string(header->numHashes));

	if (header->numBits < 2)
		fatal ("error: BloomFilter::identify_content(" + filename + ")"
		       " too few bits in vector: " + std::to_string(header->numBits));
	if (header->hashModulus < header->numBits)
		fatal ("error: BloomFilter::identify_content(" + filename + ")"
		       " hash modulus (" + std::to_string(header->hashModulus) + ")"
		     + " is less than bits in vector (" + std::to_string(header->numBits) + ")");

	// extract the info for each bitvector

	vector<BitVectorInfo> bvInfoList;

	u64 expectedOffset = header->headerSize;
	for (u32 bvIx=0 ; bvIx<header->numVectors ; bvIx++)
		{
		BitVectorInfo bvInfo;
		bvInfo.compressor = header->info[bvIx].compressor;
		bvInfo.offset     = header->info[bvIx].offset;
		bvInfo.numBytes   = header->info[bvIx].numBytes;
		u32 nameOffset    = header->info[bvIx].name;
		// we'll copy header->info[bvIx].filterInfo when we create the bit vector

		if (bvInfo.offset < header->headerSize)
			fatal ("error: BloomFilter::identify_content(" + filename + ")"
			       " offset to bitvector-" + std::to_string(bvIx)
			     + " data is within header: " + std::to_string(bvInfo.offset));

		if (bvInfo.offset != expectedOffset)
			fatal ("error: BloomFilter::identify_content(" + filename + ")"
			       " offset to bitvector-" + std::to_string(bvIx)
			     + " is " + std::to_string(bvInfo.offset)
			     + " but we expected it to be " + std::to_string(expectedOffset));

		if (nameOffset >= header->headerSize)
			fatal ("error: BloomFilter::identify_content(" + filename + ")"
			       " offset to bitvector-" + std::to_string(bvIx)
			     + " name is beyond header: " + std::to_string(nameOffset));

		u32 rrrBlockSize, rrrRankPeriod;
		switch (bvInfo.compressor & 0x000000FF)
			{
			case bvcomp_uncompressed:
			case bvcomp_roar:
			case bvcomp_unc_roar:
			case bvcomp_zeros:
			case bvcomp_ones:
				if ((bvInfo.compressor & 0xFFFFFF00) != 0) goto bad_compressor_code;
				break;
			case bvcomp_rrr:
			case bvcomp_unc_rrr:
				if ((bvInfo.compressor & 0xFF000000) != 0) goto bad_compressor_code;
				rrrBlockSize  = (bvInfo.compressor >> 8)  & 0x000000FF;
				rrrRankPeriod = (bvInfo.compressor >> 16) & 0x000000FF;
				if (rrrBlockSize != RRR_BLOCK_SIZE)
					fatal ("error: BloomFilter::identify_content(" + filename + ")"
					       " bitvector-" + std::to_string(bvIx)
					     + ", rrr block size mismatch"
					     + "\nthe file's block size is " + std::to_string(rrrBlockSize)
					     + ", program's block size is " + std::to_string(RRR_BLOCK_SIZE)
					     + "\n(see notes regarding RRR_BLOCK_SIZE in bit_vector.h)");
				if (rrrRankPeriod == 0) rrrRankPeriod = DEFAULT_RRR_RANK_PERIOD;
				if (rrrRankPeriod != RRR_RANK_PERIOD)
					fatal ("error: BloomFilter::identify_content(" + filename + ")"
					       " bitvector-" + std::to_string(bvIx)
					     + ", rrr rank period mismatch"
					     + "\nthe file's rank period is " + std::to_string(rrrRankPeriod)
					     + ", program's rank period is " + std::to_string(RRR_RANK_PERIOD)
					     + "\n(see notes regarding RRR_RANK_PERIOD in bit_vector.h)");
				bvInfo.compressor &= 0x000000FF;
				break;
			default:
			bad_compressor_code:
				fatal ("error: BloomFilter::identify_content(" + filename + ")"
				       " bitvector-" + std::to_string(bvIx)
				     + ", bad compressor code: " + std::to_string(bvInfo.compressor));
			}

		if (nameOffset != 0)
			{
			bvInfo.name = string(((char*) header) + nameOffset);
			int whichBv = bvIx % vectorsPerFilter;
			if (whichBv != 0)
				{
				string expectedName = bvInfoList[bvIx-whichBv].name;
				if (bvInfo.name != expectedName)
					fatal ("error: BloomFilter::identify_content(" + filename + ")"
					       " bitvector-" + std::to_string(bvIx)
					     + ", name is \"" + bvInfo.name + "\""
					     + " but we expected it to be \"" + expectedName + "\"");
				}
			}
		else if (numFilters == 1)
			bvInfo.name = default_filter_name(filename);
		else
			bvInfo.name = default_filter_name(filename,bvIx);

		bvInfoList.emplace_back(bvInfo);

		expectedOffset += bvInfo.numBytes;
		}

	// create the bloom filters and bit vectors; note that this will *not*
	// usually load the bit vectors
	//
	// nota bene: we assume that for, e.g. all/some filters, the bit vectors
	//   .. are in the file in the same order as they are in the filter's bvs
	//   .. array; this *should* be the same order in which they were written
	//   .. to the file

	BloomFilter* bf = nullptr;
	int infoIx = -1;
	for (const auto& bvInfo : bvInfoList)
		{
		int whichBv = (++infoIx) % vectorsPerFilter;

		if (whichBv == 0)
			{
			bf = bloom_filter(header->bfKind,
			                  filename, header->smerSize,
			                  header->numHashes, header->hashSeed1, header->hashSeed2,
			                  header->numBits, header->hashModulus);

			if (prefix.version == bffileheaderVersion1)
				{
				bf->setSizeKnown = false;
				bf->setSize      = 0;
				}
			else
				{
				if (header->setSizeKnown > 1)
					fatal ("error: BloomFilter::identify_content(" + filename + ")"
					       " set size known flag (" + std::to_string(header->setSizeKnown) + ")"
					     + " is not zero or one");
				bf->setSizeKnown = (header->setSizeKnown == 1);
				bf->setSize      = header->setSize;
				}

			bf->bvs[0] = BitVector::bit_vector(filename,bvInfo.compressor,bvInfo.offset,bvInfo.numBytes);
			}
		else
			{
			bf->bvs[whichBv] = BitVector::bit_vector(filename,bvInfo.compressor,bvInfo.offset,bvInfo.numBytes);
			}

		bf->bvs[whichBv]->filterInfo = header->info[infoIx].filterInfo;

		if (whichBv == vectorsPerFilter-1)
			{
			bf->ready = true; //…… do we really want to do this here?
			content.emplace_back(bvInfo.name,bf);
			}
		}


	if (header != nullptr) delete[] header;

	return content;
	}

//----------
//
// false_positive_rate--
//	Estimate the smer false positive rate of a bloom filter.
//
//----------
//
// Arguments:
//	u32		numHashes:	The bloom filter's number of hash functions.
//	u64		numBits:	The number of bits in the bloom filter's (conceptual)
//						bit vector.
//	u64		numItems:	The number of distinct items that have been inserted
//						into the bloom filter. This is the size of the set the
//						bloom filter represents.
//
// Returns:
//	An estimate of the false positive rate.
//
//----------
//
// Implementation notes:
//	(1)	The formula comes from the bloom filter wikipedia page (reference [2]).
//		In that description, k is the number of hash functions, m is the number
//		of bits, and n is the number of inserted elements.
//
//----------

double BloomFilter::false_positive_rate
   (u32		numHashes,
	u64		numBits,
	u64		numItems)
	{
	double fpRate;
	if (numHashes < 1)
		fpRate = 1.0;
	else if (numHashes == 1)
		fpRate = 1 - exp(-double(numItems)/double(numBits));
	else
		fpRate = pow(1 - exp(-double(numHashes*numItems)/double(numBits)),numHashes);
	return fpRate;
	}

