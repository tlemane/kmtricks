// query.cc-- classes representing queries.

#include <string>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <algorithm>

#include "utilities.h"
#include "query.h"

using std::string;
using std::vector;
using std::pair;
using std::set;
using std::cout;
using std::cerr;
using std::endl;
#define u32 std::uint32_t
#define u64 std::uint64_t

//----------
//
// Query--
//
//----------

Query::Query
   (const querydata& qd,
	double _threshold)
	  :	threshold(_threshold),
		numHashes(0),
		neededToPass(0),
		neededToFail(0),
		numUnresolved(0),
		numPassed(0),
		numFailed(0)
	{
	batchIx = qd.batchIx;
	name    = qd.name;
	seq     = qd.seq;
	}

Query::Query
   (const querydata& qd,
	double _threshold,
    std::shared_ptr<km::Repartition> rep,
    std::shared_ptr<km::HashWindow> hashwin,
    uint32_t minimsize)
	  :	threshold(_threshold),
		numHashes(0),
		neededToPass(0),
		neededToFail(0),
		numUnresolved(0),
		numPassed(0),
		numFailed(0),
        m_repartitor(rep),
        m_hash_win(hashwin),
        m_minim_size(minimsize)
	{
	batchIx = qd.batchIx;
	name    = qd.name;
	seq     = qd.seq;
	}

Query::~Query()
	{
	}


void Query::smerize
   (BloomFilter*	bf)
	{
	bf->preload(); // $$$ pierre: why preload the bf for each query ?
	u32 smerSize = bf->smerSize;

	if (bf->numHashes > 1)
		fatal ("internal error: "
		     + bf->identity() + " uses more than one hash function");

	smerHashes.clear();

	// if the sequence is too short, there are no smers

	if (seq.length() < smerSize)
		return;

	// scan the sequence's smers, convert to hash positions, and collect the
	// distinct positions; optionally collect the corresponding smers

	pair<set<u64>::iterator,bool> status;

	size_t goodNtRunLen = 0;

	for (size_t ix=0; ix<seq.length() ; ix++)
		{
		if (not nt_is_acgt(seq[ix])) { goodNtRunLen = 0;  continue; }
		if (++goodNtRunLen < smerSize) continue;

		string mer = seq.substr(ix+1-smerSize,smerSize);
        u64 hash_value = 0;
        if (m_repartitor)
        {
          km::const_loop_executor<0, KMER_N>::exec<KmerHash>(smerSize, mer, m_hash_win, m_repartitor, m_minim_size, hash_value);
        }
        else
        {
		  hash_value = bf->mer_to_hash_value(mer);
        }
		if (hash_value != BloomFilter::npos)
			{
			smerHashes.emplace_back(std::pair<std::uint64_t, std::size_t>(hash_value, ix - smerSize + 1));
			}
		}
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

void Query::read_query_file
   (std::istream&	in,
	const string&	_filename,
	double			threshold,
	vector<Query*>&	queries,
    std::string& repartFileName,
    std::string& winFileName)
	{
	bool			fileTypeKnown = false;
	bool			haveFastaHeaders = false;
	querydata		qd;

	// derive a name to use for nameless sequences

    std::shared_ptr<km::Repartition> repartitor = std::make_shared<km::Repartition>(repartFileName, "");
    std::shared_ptr<km::HashWindow> hwin = std::make_shared<km::HashWindow>(winFileName);
	string filename(_filename);
	if (filename.empty())
		filename = "(stdin)";

	string baseName(strip_file_path(_filename));

	if ((is_suffix_of (baseName, ".fa"))
	 || (is_suffix_of (baseName, ".fasta")))
		{
		string::size_type dotIx = baseName.find_last_of(".");
		baseName = baseName.substr(0,dotIx);
		}

	if (baseName.empty())
		baseName = "query";

	// read the sequences

	qd.name = "";

	string line;
	int lineNum = 0;
	int queryLineNum = 0;
	while (std::getline (in, line))
		{
		lineNum++;
		if (line.empty()) continue;

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
				fatal ("sequences precede first fasta header in \"" + filename + "\""
				     + " (at line " + std::to_string(lineNum) + ")");
			if (not qd.name.empty())
				{
				if (qd.seq.empty())
					cerr << "warning: ignoring empty sequence in \"" << filename << "\""
					     << " (at line " << std::to_string(queryLineNum) << ")" << endl;
				else
					{
					qd.batchIx = queries.size();
					queries.emplace_back(new Query(qd, threshold, repartitor, hwin, hwin->minim_size()));
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
			qd.seq += line; // $$$ pierre: why to store the sequence ?
			}

		// otherwise we're in line-by-line mode, add this line to the list

		else
			{
			qd.batchIx = queries.size();
			qd.name = baseName + std::to_string(lineNum);
			qd.seq  = line;
			queries.emplace_back(new Query(qd,threshold, repartitor, hwin, hwin->minim_size()));
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
			queries.emplace_back(new Query(qd,threshold, repartitor, hwin, hwin->minim_size()));
			}
		}

	}

