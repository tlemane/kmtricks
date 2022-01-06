#include <stdint.h>
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>
#include <utility>

#include "encoding.hpp"
#include "kff_io.hpp"

/* TO READ BEFORE ANY USAGE !!
This file provide function for dna sequences represented as a uint8_t array.
The arrays are big endian represented.
The begining of the prefix is present at the byte 0.
Byte 0 also have useless bits at the begining if the sequence size is not 4 divisable.
*/

#ifndef SEQ_H
#define SEQ_H

/** Generic object to read sequence files.
 * Each call to next_sequence should return a sequence
 **/
class SequenceStream {
public:
  /** Read the next sequence in the file.
    * If this sequence is larger than max_seq_size, it will fill seq with the max_seq_size first
    * nucleotides and return the full sequence size. remaining nucleotides are discarded.
    * @param seq Pointer that will be filled with the sequence array. The array is erased during the
    * next call to the method.
    * @return size of the sequence that have been read.
    **/
  virtual int next_sequence(uint8_t * & seq, uint max_seq_size, uint8_t * & data, uint max_data_size) = 0;
};


/** Read txt sequence file (1 sequence per line)
 */
class KffSeqStream : public SequenceStream {
private:
public:
  Kff_reader reader;
  KffSeqStream(const std::string filename)
      : reader(filename)
  {};
  /** Read the next sequence (the whole block) storing the nucleotides in the seq array and the
   * associated data in the data array. These two arrays must be pre-allocated.
   * 
   * @param seq User allocated array where the sequence will be loaded
   * @param max_seq_size Maximum sequence size that can fit in the seq buffer (in nucleotides).
   * @param data User allocated array where the data will be loaded
   * @param max_data_size Maximum number of Bytes that can fit into data.
   * 
   * @return The number of kmers into the sequence. -1 if the sequence do not fit in seq or the
   * data in data.
   **/
  int next_sequence(uint8_t * & seq, uint max_seq_size, uint8_t * & data, uint max_data_size);
};


/** Extract a subsequence of sequence.
  * @param sequence Original sequence
  * @param seq_size Size in nucleotides of the sequence
  * @param extracted A memory space already allocated by the user to copy the subsequence.
  * This memory must be 1 byte larger than the requiered space for the subsequence.
  * @param begin_nucl first nucleotide to extract (between 0 and seq_size-1)
  * @param end_nucl last nucleotide to extract (between 0 and seq_size-1)
  */
void subsequence(const uint8_t * sequence, const uint seq_size, uint8_t * extracted, const uint begin_nucl, const uint end_nucl);


/** Compare two subsequences. -1 if the first one is smaller in alpha order +1 is the second one
  * 0 if equals
  * 
  * 
  */
int sequence_compare(const uint8_t * seq1, const uint seq1_size,
                      const uint seq1_start, const uint seq1_stop,
                      const uint8_t * seq2, const uint seq2_size,
                      const uint seq2_start, const uint seq2_stop);


/** Translate a sequence to an 64 bits integer. If the sequence length is more than 32, then only
  * the last 32 nucleotides are used for the conversion.
  *
  * @param seq A nucleotide sequence to translate
  * @param seq_size Size in nucleotides of the sequence
  */
uint64_t seq_to_uint(const uint8_t * seq, uint seq_size);

/** Translate a subsequence into a 64 bits integer. If the sequence size > 32, then the 32 suffix
  * of the sequence is used.
  * 
  * @param seq The sequence sur translate
  * @param seq_size The number of nucleotides in the sequence
  * @param start_nucl First nucletide index of the target subsequence
  * @param end_nucl Last nucleotide of the subsequence to translate
  * 
  * @return Translate subsequence
  */
uint64_t subseq_to_uint(const uint8_t * seq, uint seq_size, uint start_nucl, uint end_nucl);

/** Translate a binarized uint sequence into an binarized sequence array.
  * @param seq Sequence stored in a uint (ie max 32 nucleotides)
  * @param bin_seq A Byte array to store the sequence. Must be allocated.
  * @param size Number of nucleotides in the sequence.
  */
void uint_to_seq(uint seq, uint8_t * bin_seq, uint size);


// ----- Minimizer search related functions -----

typedef struct {
  uint64_t start_position;
  uint64_t stop_position;
  int64_t minimizer_position;
  uint64_t minimizer;
} skmer;

class MinimizerSearcher {
public:
  uint k;
  uint m;
  uint add_count;
  uint use_count;
  uint max_seq_size;
  bool single_side;
  std::vector<uint64_t> candidates;
  std::vector<bool> is_rev_candidates;
  std::vector<uint64_t> mini_buffer;
  std::vector<uint64_t> minis;
  std::vector<int64_t> mini_pos;
  std::vector<std::pair<uint64_t, uint64_t> > skmers;
  uint64_t nucl_fwd[4][256];
  uint64_t nucl_rev[4][256];
  RevComp rc;
  MinimizerSearcher(const uint k, const uint m, const uint8_t encoding[4], const uint max_seq_size = 0, const bool single_side = false)
          : k(k), m(m), add_count(0), use_count(0), max_seq_size(max_seq_size), single_side(single_side)
          , mini_buffer(max_seq_size < m - 1 ? 0 : (max_seq_size - m + 1) * 2, 0)
          , minis(max_seq_size < k - 1 ? 0 : max_seq_size - k + 1, 0)
          , mini_pos(max_seq_size < k - 1 ? 0 : max_seq_size - k + 1, 0)
          , skmers()
          , rc(encoding)
  {
    for (uint byte=0 ; byte<256 ; byte++) {
      for (uint nucl_pos=0 ; nucl_pos<4 ; nucl_pos++) {
        nucl_fwd[nucl_pos][byte] = (byte >> (2 * (3 - nucl_pos))) & 0b11;
        nucl_rev[nucl_pos][byte] = this->rc.reverse[nucl_fwd[nucl_pos][byte]] << (2 * (this->m - 1));
      }
    }
  };

  ~MinimizerSearcher() {
    // std::cout << add_count << std::endl << use_count << std::endl << std::endl;
  }

  /** Fill the first half of the mini_buffer with m-mers candidates for the fwd.
   * Same with the second half and the candidates from the rev-comp.
   * 
   * @param seq Binarized sequence
   * @param seq_size Sequence size
   **/
  void compute_candidates(const uint8_t * seq, const uint seq_size);
  /** Fill the mini_pos vector with one sequence position per kmer.
   * Positive/zero number mean on forward, negative on reverse (complement a 1 to remove ambiguity of +0 and -0).
   * 1 means that the minimizer starts at position 1 on the forward sequence.
   * -3 means that the mini starts at position 2 (comp a 1) on the sequence and have to be reverse complemented.
   * 
   * @param nb_kmers The number of kmers inside the sequence
   **/
  void compute_minimizers(const uint nb_kmers);
  /** Compute the boundaries of the superkmer on the fwd strand
   * 
   * @param nb_kmers Number of kmers in the sequence
   **/
  void compute_skmers(const uint nb_kmers);

  /** Get a vector of all the skmers in a sequence.
   * All the other methods of this class are called by this function.
   * No precomputation needed.
   * 
   * @param seq Sequence to analyse
   * @param seq_size Sequence size in nucleotides.
   * 
   * @return A vector containing object of type skmer.
   **/
  std::vector<skmer> get_skmers(const uint8_t * seq, const uint seq_size);
  std::vector<skmer> get_skmers_fast(const uint8_t * seq, const uint seq_size);
};




// ----- Usefull binary functions -----
void leftshift8(uint8_t * bitarray, size_t length, size_t bitshift);
void rightshift8(uint8_t * bitarray, size_t length, size_t bitshift);
uint8_t fusion8(uint8_t left_bits, uint8_t right_bits, size_t merge_index);

#endif