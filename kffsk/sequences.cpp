#include <cstring>
#include <cassert>
#include <algorithm>

#include "sequences.hpp"



using namespace std;


/* Bitshift to the left all the bits in the array with a maximum of 7 bits.
 * Overflow on the left will be set into the previous cell.
 */
void leftshift8(uint8_t * bitarray, size_t length, size_t bitshift) {
	assert(bitshift < 8);

	for (uint64_t i=0 ; i<length-1 ; i++) {
		bitarray[i] = (bitarray[i] << bitshift) | (bitarray[i+1] >> (8-bitshift));
	}
	bitarray[length-1] <<= bitshift;
}

/* Similar to the previous function but on the right */
void rightshift8(uint8_t * bitarray, size_t length, size_t bitshift) {
	assert(bitshift < 8);

	for (uint64_t i=length-1 ; i>0 ; i--) {
		bitarray[i] = (bitarray[i-1] << (8-bitshift)) | (bitarray[i] >> bitshift);
	}
	bitarray[0] >>= bitshift;
}

/* Fusion to bytes into one.
 * The merge_index higher bits are from left_bits the others from right_bits
 */
uint8_t fusion8(uint8_t left_bits, uint8_t right_bits, size_t merge_index) {
	uint8_t mask = 0xFF << (8-merge_index);
	return (left_bits & mask) | (right_bits & ~mask);
}

int KffSeqStream::next_sequence(uint8_t * & seq, uint max_seq_size, uint8_t * & data, uint max_data_size) {
// int KffSeqStream::next_sequence(uint8_t * & seq, uint max_seq_size, uint8_t * & data, uint max_data_size) {
	if (this->reader.has_next()) {
		uint max_seq = this->reader.k + this->reader.max - 1;
		uint max_data = this->reader.max * this->reader.data_size;
		if (max_seq_size < max_seq or max_data_size < max_data)
			return -1;
		else
			return this->reader.next_block(seq, data);
	}

	return 0;
}




// #include <iostream>
void subsequence(const uint8_t * sequence, const uint seq_size, uint8_t * extracted, const uint begin_nucl, const uint end_nucl) {
	// Extract the correct slice
	uint seq_left_offset = (4 - seq_size % 4) % 4;
	uint extract_start_byte = (seq_left_offset + begin_nucl) / 4;
	uint extract_stop_byte = (seq_left_offset + end_nucl) / 4;

	memcpy(extracted, sequence + extract_start_byte, extract_stop_byte - extract_start_byte + 1);

	// Align the bits
	uint extract_left_offset = (seq_left_offset + begin_nucl) % 4;
	uint extract_right_offset = (seq_size - end_nucl - 1) % 4;
	

	if (extract_right_offset < 4 - extract_left_offset) {
		rightshift8(extracted, extract_stop_byte - extract_start_byte + 1, extract_right_offset * 2);
	} else {
		leftshift8(extracted, extract_stop_byte - extract_start_byte + 1, (4 - extract_right_offset) * 2);
	}
}


int sequence_compare(const uint8_t * seq1, const uint seq1_size,
											const uint seq1_start, const uint seq1_stop,
											const uint8_t * seq2, const uint seq2_size,
											const uint seq2_start, const uint seq2_stop) {
	// If inequal size
	if (seq1_stop - seq1_start != seq2_stop - seq2_start)
		return (seq1_stop - seq1_start) < (seq2_stop - seq2_start) ? -1 : 1;

	// Extraction of subsequences
	uint subseq_size = seq1_stop - seq1_start + 1;
	uint subseq_bytes = (subseq_size + 3) / 4;
	uint8_t * subseq1 = new uint8_t[subseq_bytes];
	subsequence(seq1, seq1_size, subseq1, seq1_start, seq1_stop);
	uint8_t * subseq2 = new uint8_t[subseq_bytes];
	subsequence(seq2, seq2_size, subseq2, seq2_start, seq2_stop);

	// comparison of the subsequences (same size)
	uint return_val = 0;
	// First byte
	uint offset = (4 - (subseq_size % 4)) % 4;
	uint mask = (1 << (2 * (4 - offset))) - 1;
	if ((subseq1[0] & mask) != (subseq2[0] & mask))
		return_val = (subseq1[0] & mask) < (subseq2[0] & mask) ? -1 : 1;

	// All following bytes
	for (uint b=1 ; b<subseq_bytes and return_val==0 ; b++) {
		if (subseq1[b] != subseq2[b]) {
			return_val = subseq1[b] < subseq2[b] ? -1 : 1;
		}
	}

	delete[] subseq1;
	delete[] subseq2;
	return return_val;
}


void MinimizerSearcher::compute_candidates(const uint8_t * seq, const uint seq_size) {
	if (seq_size > this->max_seq_size) {
		this->max_seq_size = seq_size;
		this->mini_buffer.resize((max_seq_size - m + 1) * 2, 0);
		this->minis.resize(this->max_seq_size - k + 1, 0);
		this->mini_pos.resize(this->max_seq_size - k + 1, 0);
	}

	uint offset = (4 - (seq_size % 4)) % 4;
	uint64_t current_value = 0;
	uint64_t current_rev_value = 0;
	// Compute prefix of first candidate
	for (uint i=0 ; i<this->m-1 ; i++) {
		uint idx = offset + i;
		uint byte_idx = idx/4;
		uint nucl_shift = 3 - (idx % 4);

		uint nucl = (seq[byte_idx] >> (nucl_shift * 2)) & 0b11;
		current_value = (current_value << 2) + nucl;
		current_rev_value = (current_rev_value >> 2) + (this->rc.reverse[nucl] << (2 * (this->m - 1)));
	}
	
	// Compute minimizer candidates
	uint64_t m_mask = (1 << (this->m*2)) - 1;
	for (uint i=this->m-1, kmer_idx=0 ; i<seq_size ; i++, kmer_idx++) {
		uint idx = offset + i;
		uint byte_idx = idx/4;
		uint nucl_shift = 3 - (idx % 4);

		uint nucl = (seq[byte_idx] >> (nucl_shift * 2)) & 0b11;
		current_value = ((current_value << 2) + nucl) & m_mask;
		current_rev_value = (current_rev_value >> 2) + this->nucl_rev[idx%4][seq[byte_idx]];//(this->rc.reverse[nucl] << (2 * (m - 1)));
		this->mini_buffer[kmer_idx] = current_value;
		this->mini_buffer[this->mini_buffer.size()/2 + kmer_idx] = current_rev_value;
	}
}



// void MinimizerSearcher::compute_minimizers(const uint nb_kmers) {
// 	// Compute the minimizer of each sliding window of size k - m
// 	uint max_nb_candidates = this->mini_buffer.size()/2;
// 	for (uint i=0 ; i<max_nb_candidates ; i++) {
// 		if (i < this->minis.size())
// 			this->minis[i] = 0xFFFFFFFFFFFFFFFF;
// 		// Compute minimum
// 		uint64_t mini = this->mini_buffer[i];
// 		int64_t abs_pos = i;
// 		if (not this->single_side) {
// 			if (this->mini_buffer[max_nb_candidates + i] < mini) {
// 				mini = this->mini_buffer[max_nb_candidates + i];
// 				abs_pos = - (int)i - 1;
// 			}
// 		}

// 		// Fill the minis
// 		for (int j=max(0, (int)i-(int)k+(int)m) ; j<(int)this->minis.size() ; j++) {
// 			if (mini < this->minis[j]) {
// 				this->minis[j] = mini;
// 				this->mini_pos[j] = abs_pos;
// 			}
// 		}
// 	}
// }

void MinimizerSearcher::compute_minimizers(const uint nb_kmers) {
	// Compute the minimizer of each sliding window of size k - m
	uint max_nb_candidates = this->mini_buffer.size()/2;
	for (uint i=0 ; i<nb_kmers ; i++) {
		auto smallest_fwd = min_element(
			this->mini_buffer.begin() + i,
			this->mini_buffer.begin()+i+(k-m)+1
		);
		auto smallest_rev = smallest_fwd;
		if (not this->single_side) {
			smallest_rev = min_element(
				this->mini_buffer.begin()+max_nb_candidates+i,
				this->mini_buffer.begin()+max_nb_candidates+i+(k-m)+1
			);
		}

		int mini_pos;
		if (this->single_side or (*smallest_fwd) < (*smallest_rev)) {
			mini_pos = smallest_fwd - this->mini_buffer.begin();
		} 
		else if ((*smallest_fwd) == (*smallest_rev)) {
			if (smallest_fwd - smallest_rev <= 0) {
				mini_pos = smallest_fwd - this->mini_buffer.begin();
			} else {
				mini_pos = - (smallest_rev - (this->mini_buffer.begin() + this->mini_buffer.size() / 2)) - 1;
			}
		}
		else {
			mini_pos = - (smallest_rev - (this->mini_buffer.begin() + this->mini_buffer.size() / 2)) - 1;
		}

		this->mini_pos[i] = mini_pos;
	}
}

void MinimizerSearcher::compute_skmers(const uint nb_kmers) {
	this->skmers.clear();

	uint last_mini_start = 0;
	for (uint idx=1 ; idx<nb_kmers ; idx++) {
		if (this->mini_pos[idx] != this->mini_pos[last_mini_start]) {
			this->skmers.emplace_back(last_mini_start, idx - 1 + this->k - 1);
			last_mini_start = idx;
		}
	}
	this->skmers.emplace_back(last_mini_start, nb_kmers - 1 + this->k - 1);
}


vector<skmer> MinimizerSearcher::get_skmers(const uint8_t * seq, const uint seq_size) {
	// this->get_skmers_fast(seq, seq_size);

	this->compute_candidates(seq, seq_size);
	uint nb_kmers = seq_size - k + 1;
	this->compute_minimizers(nb_kmers);
	this->compute_skmers(nb_kmers);

	vector<skmer> skmers(this->skmers.size(), {0,0,0,0});

	for (uint sk_idx=0 ; sk_idx<this->skmers.size() ; sk_idx++) {
		uint64_t start = this->skmers[sk_idx].first;
		skmers[sk_idx].start_position = start;
		skmers[sk_idx].stop_position = this->skmers[sk_idx].second;
		int64_t mini_pos = this->mini_pos[start];
		skmers[sk_idx].minimizer_position = mini_pos;
		if (mini_pos >= 0)
			skmers[sk_idx].minimizer = this->mini_buffer[mini_pos];
		else {
			mini_pos = - mini_pos - 1;
			skmers[sk_idx].minimizer = this->mini_buffer[this->mini_buffer.size() / 2 + mini_pos];
		}
	}

	return skmers;
}


vector<skmer> MinimizerSearcher::get_skmers_fast(const uint8_t * seq, const uint seq_size) {
	vector<skmer> skmers;

	if (seq_size > this->candidates.size()) {
		this->candidates.resize(seq_size, 0);
		this->is_rev_candidates.resize(seq_size, 0);
	}

	// uint8_t encoding[] = {0, 1, 3, 2};
	// Stringifyer strif(encoding);
	// cout << strif.translate(seq, seq_size) << endl;

	// Usefull variables
	uint idx_offset = (4 - (seq_size % 4)) % 4;
	uint8_t current_byte = seq[0];
	uint64_t mini_mask = (1 << (this->m * 2)) - 1;

	// --- Preprocess the m-1 size first minimizer prefix ---
	uint64_t current_candidate = 0;
	uint64_t current_candidate_fwd = 0;
	uint64_t current_candidate_rev = 0;
	bool current_rev = false;
	for (uint seq_idx=0 ; seq_idx<this->m-1 ; seq_idx++) {
		uint64_t abs_idx = idx_offset + seq_idx;
		// Load new byte
		uint64_t nucl_idx = abs_idx % 4;
		if (nucl_idx == 0)
			current_byte = seq[abs_idx / 4];

		// Get nucleotide
		// uint64_t nucl = (current_byte >> ((3 - nucl_idx) * 2)) & 0b11;
		// Grow the minimizers
		current_candidate_fwd <<= 2;
		current_candidate_fwd |= this->nucl_fwd[nucl_idx][current_byte];
		if (not this->single_side) {
			current_candidate_rev >>= 2;
			current_candidate_rev |= this->nucl_rev[nucl_idx][current_byte];
		}
		// cout << strif.translate(current_candidate, seq_idx+1) << endl;
	}

	uint64_t current_minimizer = 0xFFFFFFFFFFFFFFFF;
	uint64_t abs_mini_pos = 0;
	bool mini_rev = false;

	// --- Process the minimizer for the first kmer ---
	for (uint seq_idx = this->m-1 ; seq_idx<this->k ; seq_idx++) {
		// --- Current candidate computation ---
		uint64_t abs_idx = idx_offset + seq_idx;
		// Load new byte
		uint64_t nucl_idx = abs_idx % 4;
		if (nucl_idx == 0)
			current_byte = seq[abs_idx / 4];

		// Get nucleotide
		// uint64_t nucl = (current_byte >> ((3 - nucl_idx) * 2)) & 0b11;
		// Update current minimizer fwd
		current_candidate_fwd <<= 2;
		current_candidate_fwd |= this->nucl_fwd[nucl_idx][current_byte];
		current_candidate_fwd &= mini_mask;
		current_candidate = current_candidate_fwd;

		if (not this->single_side) {
			current_candidate_rev >>= 2;
			current_candidate_rev |= this->nucl_rev[nucl_idx][current_byte];
			current_rev = false;
			if (current_candidate_rev < current_candidate_fwd) {
				current_candidate = current_candidate_rev;
				current_rev = true;
			}
		}

		// current minimizer init
		if (current_candidate < current_minimizer) {
			current_minimizer = current_candidate;
			abs_mini_pos = seq_idx - (this->m - 1);
			mini_rev = current_rev;
		}
	}

	skmer sk = {0, 0, 0, 0};

	// --- Compute all the other kmers ---
	for (uint seq_idx = this->k ; seq_idx<seq_size ; seq_idx++) {
		// --- Current candidate computation ---
		uint64_t abs_idx = idx_offset + seq_idx;
		// Load new byte
		uint64_t nucl_idx = abs_idx % 4;
		if (nucl_idx == 0)
			current_byte = seq[abs_idx / 4];

		// Get nucleotide
		// uint64_t nucl = (current_byte >> ((3 - nucl_idx) * 2)) & 0b11;
		// Update current minimizer fwd
		current_candidate_fwd <<= 2;
		current_candidate_fwd |= this->nucl_fwd[nucl_idx][current_byte];
		current_candidate_fwd &= mini_mask;
		current_candidate = current_candidate_fwd;
		// Update current mini rev
		if (not this->single_side) {
			current_candidate_rev >>= 2;
			current_candidate_rev |= this->nucl_rev[nucl_idx][current_byte];;
			current_rev = false;
			if (current_candidate_rev < current_candidate_fwd) {
				current_candidate = current_candidate_rev;
				current_rev = true;
			}
		}

		bool mini_change = false;

		// cout << seq_idx - m + 1 << " " << seq_idx - k + 1 << "\t" << strif.translate(current_candidate, this->m) << " " << current_candidate;
		// --- get the best minimizer for current kmer ---
		// New minimizer found
		if (current_candidate < current_minimizer) {
			sk.minimizer = current_minimizer;
			sk.minimizer_position = abs_mini_pos;
			if (mini_rev)
				sk.minimizer_position = - sk.minimizer_position - 1;

			current_minimizer = current_candidate;
			abs_mini_pos = seq_idx - (this->m - 1);
			mini_rev = current_rev;
			// cout << " X (" << abs_mini_pos << ")";
			mini_change = true;
		}
		// Previous mini OUTDATED
		else if (seq_idx == abs_mini_pos + k) {
			cerr << "Not properly implemented. Quitting" << endl;
			exit(1);
			sk.minimizer = current_minimizer;
			sk.minimizer_position = abs_mini_pos;
			if (mini_rev)
				sk.minimizer_position = - sk.minimizer_position - 1;

			// cout << " OUTDATED ";
			current_minimizer = candidates[(seq_idx-k+m)];
			this->use_count += 1;
			abs_mini_pos = seq_idx - k + 1;
			mini_rev = current_rev = is_rev_candidates[(seq_idx-k+m)];
			// Recompute from the saved candidates
			for (uint idx=seq_idx-k+m+1 ; idx<=seq_idx ; idx++) {
				this->use_count += 1;
				if (candidates[idx] < current_minimizer) {
					current_minimizer = candidates[idx];
					abs_mini_pos = idx - (this->m - 1);
					mini_rev = current_rev = is_rev_candidates[idx];
				}
			}

			mini_change = true;
			// cout << abs_mini_pos << " ";

			// cout << strif.translate(current_minimizer, this->m) << " (" << abs_mini_pos << ")";
		}
		// Keep the same minimizer
		else {}

		// --- Skmer generation ---
		if (mini_change) {
			sk.stop_position = seq_idx - 1;
			skmers.push_back(sk);
			sk.start_position = seq_idx - this->k + 1;
		}
	}

	// Last skmer
	sk.stop_position = seq_size - 1;
	sk.minimizer = current_minimizer;
	sk.minimizer_position = abs_mini_pos;
	if (mini_rev)
		sk.minimizer_position = - sk.minimizer_position - 1;

	skmers.push_back(sk);

	// uint8_t subseq[10];
	// subsequence(seq, seq_size, subseq, skmer_start, skmer_stop);
	// cout << " " << strif.translate(subseq, skmer_stop - skmer_start + 1) << endl << endl;

	// for (const skmer sk : skmers) {
	// 	uint8_t subseq[10];
	// 	subsequence(seq, seq_size, subseq, sk.start_position, sk.stop_position);
	// 	cout << strif.translate(subseq, sk.stop_position - sk.start_position + 1) << endl;
	// }

	// exit(0);
	return skmers;
}


uint64_t seq_to_uint(const uint8_t * seq, uint seq_size) {
	uint nucl_to_extract = seq_size;
	if (nucl_to_extract > 32)
		nucl_to_extract = 32;

	uint seq_offset = (4 - (seq_size % 4)) % 4;
	uint seq_bytes = (seq_size + 3) / 4;
	uint useless_seq_nucl = seq_size - nucl_to_extract;

	uint suff_offset = (4 - (nucl_to_extract % 4)) % 4;
	uint mask = (1 << (2 * (4 - suff_offset))) - 1;
	uint suff_first_byte = (seq_offset + useless_seq_nucl) / 4;

	uint64_t val = seq[suff_first_byte] & mask;
	for (uint idx=suff_first_byte+1 ; idx<seq_bytes ; idx++) {
		val <<= 8;
		val += seq[idx];
	}

	return val;
}


uint64_t subseq_to_uint(const uint8_t * seq, uint seq_size, uint start_nucl, uint end_nucl) {
	// Trunkate too long sequences
	if (end_nucl - start_nucl + 1 > 32)
		start_nucl = end_nucl - 31;

	// Determine boundary bytes
	uint seq_offset = (4 - (seq_size % 4)) % 4;
	uint first_sub_byte = (seq_offset + start_nucl) / 4;
	uint last_sub_byte = (seq_offset + end_nucl) / 4;

	if (first_sub_byte == last_sub_byte) {
		uint64_t val = seq[first_sub_byte];
		uint shift = 2 * (3 - seq_offset - end_nucl);
		uint mask = (1 << ((end_nucl - start_nucl + 1) * 2)) - 1;
		return (val >> shift) & mask;
	}

	// First byte integration
	uint mask = (1 << (2 * (4 - ((seq_offset + start_nucl) % 4)))) - 1;
	uint64_t sub_val = seq[first_sub_byte] & mask;

	// Middle bytes
	for (uint b=first_sub_byte+1 ; b<last_sub_byte ; b++) {
		sub_val <<= 8;
		sub_val += seq[b];
	}

	// End byte
	uint last_shift = (seq_size - end_nucl - 1) % 4;
	uint8_t end_byte = seq[last_sub_byte] >> (2 * last_shift);
	sub_val <<= 2 * (4 - last_shift);
	sub_val += end_byte;

	return sub_val;
}


void uint_to_seq(uint seq, uint8_t * bin_seq, uint size) {
	uint seq_bytes = (size + 3) / 4;

	for (int idx=seq_bytes-1 ; idx>=0 ; idx--) {
		bin_seq[idx] = seq & 0b11111111;
		seq >>= 8;
	}
}