// query.cc-- classes representing queries.

#include <string>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "utilities.h"
#include "query.h"

using std::cerr;
using std::cout;
using std::endl;
using std::pair;
using std::set;
using std::string;
using std::vector;
#define u32 std::uint32_t
#define u64 std::uint64_t

//----------
//
// Query--
//
//----------

Query::Query(const querydata &qd,
			 double _threshold)
	: threshold(_threshold),
	  numPositions(0),
	  neededToPass(0),
	  neededToFail(0),
	  numUnresolved(0),
	  numPassed(0),
	  numFailed(0),
	  nodesExamined(0)
{
	batchIx = qd.batchIx;
	name = qd.name;
	seq = qd.seq;
}

Query::Query(const querydata &qd,
			 double _threshold,
			 std::shared_ptr<km::Repartition> rep,
			 std::shared_ptr<km::HashWindow> hashwin,
			 uint32_t minimsize)
	: threshold(_threshold),
	  numPositions(0),
	  neededToPass(0),
	  neededToFail(0),
	  numUnresolved(0),
	  numPassed(0),
	  numFailed(0),
	  nodesExamined(0),
	  _repartitor(rep),
	  _hash_win(hashwin),
	  _minimsize(minimsize)
{
	batchIx = qd.batchIx;
	name = qd.name;
	seq = qd.seq;
}

Query::~Query()
{
}

template<size_t KSIZE>
struct KmerHash
{
	void operator()(const std::string& mer, std::shared_ptr<km::HashWindow> hw, std::shared_ptr<km::Repartition> repart, uint32_t minim, uint64_t& pos)
	{
		km::Kmer<KSIZE> kmer(mer); km::Kmer<KSIZE> cano = kmer.canonical();
		uint32_t part = repart->get_partition(cano.minimizer(minim).value());
		pos = km::KmerHashers<1>::WinHasher<KSIZE>(part, hw->get_window_size_bits())(cano);
	}
};

void Query::kmerize(BloomFilter *bf,
					bool distinct)
{
	bf->preload();
	u32 kmerSize = bf->kmerSize;

	if (bf->numHashes > 1)
		fatal("internal error: " + bf->identity() + " uses more than one hash function");

	kmerPositions.clear();
	kmerized2endpos.clear();

	// if the sequence is too short, there are no kmers

	if (seq.length() < kmerSize)
		return;

	// scan the sequence's kmers, convert to hash positions, and collect the
	// distinct positions; optionally collect the corresponding kmers

	set<u64> positionSet;
	pair<set<u64>::iterator, bool> status;

	size_t goodNtRunLen = 0;
	u64 part = 0;
	u64 pos = 0;
	u64 wsize = 0;
	u64 bval = 0;
	for (size_t ix = 0; ix < seq.length(); ix++)
	{

		// finds the first position containing a good kmer.
		if (not nt_is_acgt(seq[ix]))
		{
			goodNtRunLen = 0;
			continue;
		}
		if (++goodNtRunLen < kmerSize)
			continue;

		// ix ends a valid kmer/
		string mer = seq.substr(ix + 1 - kmerSize, kmerSize);

	  if (_repartitor)
		{
			km::const_loop_executor<0, KMER_N>::exec<KmerHash>(
				kmerSize, mer, _hash_win, _repartitor, _minimsize, pos);
			//km::Kmer<32> kmer(mer); km::Kmer<32> cano = kmer.canonical();
			//std::cout << kmer.to_string() << " " << cano.to_string() << std::endl;
			//part = _repartitor->get_partition(cano.minimizer(_minimsize).value());
			//std::cout << std::to_string(part) << std::endl;
			//pos = km::KmerHashers<1>::WinHasher<32>(part, _hash_win->get_window_size_bits())(cano);
			//std::cout << std::to_string(pos) << std::endl;
		}
		else
		{
			pos = bf->mer_to_position(mer);
		}

		if (pos != BloomFilter::npos)
		{
			if (distinct)
			{
				status = positionSet.insert(pos);
				if (status.second == false) // pos was already in the set
					continue;
			}
			kmerPositions.emplace_back(pos);
			kmerized2endpos.emplace_back(ix);
		}

	}
	// not needed anymore
	seq_length = seq.length();
	seq.clear(); 
	seq.shrink_to_fit(); 
}






//----------
//
// read_query_file--
//	Read queries from a file (names and nucleotide sequences), collecting them
//	in a list.
//
//----------
//
// Arguments:
//	istream&		in:			The file/stream to read from.
//	const string&	filename:	The name of the file. This is only used for
//								.. assigning names to unnamed queries, and for
//								.. error reports.
//	double			threshold:	'Hit' threshold for queries in this file.
//								0 < threshold <= 1
//	vector<Query*>&	queries:	List to copy queries to. Queries are appended
//								.. to this, so any contents it has prior to the
//								.. call are preserved. Note that the caller is
//								.. responsible for (eventually) deleting the
//								.. objects we add to this list.
//
// Returns:
//	(nothing)
//
//----------
//
// Notes:
//	(1)	We accept two sequence file formats. One format is fasta, which has
//		header lines beginning with '>'; fasta sequences can be broken into
//		multiple lines. The other format is one sequence per line.
//	(2)	If sequence names aren't available, we create them by appending the
//		file name's core with the line number.
//
//----------


void Query::read_query_file_km(std::istream &in,
							   const string &_filename,
							   double threshold,
							   vector<Query *> &queries,
							   string &repartFileName,
							   string &winFileName)
{
	bool fileTypeKnown = false;
	bool haveFastaHeaders = false;
	querydata qd;

	//km::Repartition* repartitor = new km::Repartition(repartFileName, "");
	std::shared_ptr<km::Repartition> repartitor = std::make_shared<km::Repartition>(repartFileName, "");
	std::shared_ptr<km::HashWindow> hwin = std::make_shared<km::HashWindow>(winFileName);

	uint32_t minimizer_size = hwin->minim_size();
	uint64_t h0, h1;
	uint32_t nb_parts;
	// derive a name to use for nameless sequences

	string filename(_filename);
	if (filename.empty())
		filename = "(stdin)";

	string baseName(strip_file_path(_filename));

	if ((is_suffix_of(baseName, ".fa")) || (is_suffix_of(baseName, ".fasta")))
	{
		string::size_type dotIx = baseName.find_last_of(".");
		baseName = baseName.substr(0, dotIx);
	}

	if (baseName.empty())
		baseName = "query";

	// read the sequences

	qd.name = "";

	string line;
	int lineNum = 0;
	int queryLineNum = 0;
	while (std::getline(in, line))
	{
		lineNum++;
		if (line.empty())
			continue;

		if (not fileTypeKnown)
		{
			haveFastaHeaders = (line[0] == '>');
			fileTypeKnown = true;
		}

		// if this is a fasta header, add the previous sequence to the list
		// and start a new one

		if (line[0] == '>')
		{
			if (not haveFastaHeaders)
				fatal("sequences precede first fasta header in \"" + filename + "\"" + " (at line " + std::to_string(lineNum) + ")");
			if (not qd.name.empty())
			{
				if (qd.seq.empty())
					cerr << "warning: ignoring empty sequence in \"" << filename << "\""
						 << " (at line " << std::to_string(queryLineNum) << ")" << endl;
				else
				{
					qd.batchIx = queries.size();
					queries.emplace_back(new Query(qd, threshold, repartitor, hwin, minimizer_size));
				}
			}

			queryLineNum = lineNum;
			qd.name = strip_blank_ends(line.substr(1));
			if (qd.name.empty())
				qd.name = baseName + std::to_string(lineNum);
			qd.seq = "";
		}

		// if it's not a fasta header, and we're in fasta mode, add this line
		// to the current sequence

		else if (haveFastaHeaders)
		{
			qd.seq += line;
		}

		// otherwise we're in line-by-line mode, add this line to the list

		else
		{
			qd.batchIx = queries.size();
			qd.name = baseName + std::to_string(lineNum);
			qd.seq = line;
			queries.emplace_back(new Query(qd, threshold, repartitor, hwin, minimizer_size));
			qd.name = "";
		}
	}

	// if we were accumulating a sequence, add it to the list

	if (not qd.name.empty())
	{
		if (qd.seq.empty())
			cerr << "warning: ignoring empty sequence in \"" << filename << "\""
				 << " (preceding line " << lineNum << ")" << endl;
		else
		{
			qd.batchIx = queries.size();
			queries.emplace_back(new Query(qd, threshold, repartitor, hwin, minimizer_size));
		}
	}
}
