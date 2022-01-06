#include "encoding.hpp"
#include "sequences.hpp"

using namespace std;



Translator::Translator(uint8_t source[4], uint8_t destination[4]) {
	// Nucleotide translation values
	uint8_t nucl_translation[4];
	for (uint i=0 ; i<4 ; i++)
		nucl_translation[source[i]] = destination[i] & 0b11;
	
	// Lookup translation for Bytes
	for (uint i=0 ; i<256 ; i++) {
		lookup[i] = 0;

		for (uint pos=0 ; pos<4 ; pos++) {
			// Get nucleotide
			uint8_t letter = (i >> (2*pos)) & 0b11;
			// Change encoding
			letter = nucl_translation[letter];
			// Writ in the lookup table
			lookup[i] = lookup[i] | (letter << (2*pos));
		}
	}
}

void Translator::translate(uint8_t * sequence, size_t byte_length) {
	// Translate Byte per Byte
	for (size_t idx=0 ; idx<byte_length ; idx++) {
		// Translate a Byte
		sequence[idx] = lookup[sequence[idx]];
	}
}


RevComp::RevComp(const uint8_t encoding[4]) {
	this->reverse[encoding[0]] = encoding[3];
	this->reverse[encoding[1]] = encoding[2];
	this->reverse[encoding[2]] = encoding[1];
	this->reverse[encoding[3]] = encoding[0];

	for (uint i=0 ; i<256 ; i++) {
		uint8_t val = i;
		uint8_t rc_val = 0;

		for (uint j=0 ; j<4 ; j++) {
			rc_val <<= 2;
			rc_val += this->reverse[val & 0b11];
			val >>= 2;
		}

		this->translations[i] = rc_val;
	}
}

void RevComp::rev_comp(uint8_t * seq, const uint64_t seq_size) const {
	uint nb_bytes = (seq_size + 3) / 4;
	// reverse and translate each byte
	for (uint byte_idx=0 ; byte_idx<(nb_bytes+1)/2 ; byte_idx++) {
		uint8_t save = this->translations[seq[byte_idx]];
		seq[byte_idx] = this->translations[seq[nb_bytes-1-byte_idx]];
		seq[nb_bytes-1-byte_idx] = save;
	}

	uint8_t offset = (4 - (seq_size % 4)) % 4;
	rightshift8(seq, nb_bytes, offset*2);
}


void RevComp::rev_data(uint8_t * data, const uint64_t data_size, const uint64_t nb_kmers) const {
	for (uint idx=0 ; idx<nb_kmers/2 ; idx++) {
		uint rev_idx = nb_kmers - 1 - idx;
		for (uint data_idx=0 ; data_idx<data_size ; data_idx++) {
			uint8_t save = data[idx * data_size + data_idx];
			data[idx * data_size + data_idx] = data[rev_idx * data_size + data_idx];
			data[rev_idx * data_size + data_idx] = save;
		}
	}
}


Stringifyer::Stringifyer(uint8_t encoding[4]) {
	// Nucleotide translation values
	string nucl_translation[4];
	nucl_translation[encoding[0]] = "A";
	nucl_translation[encoding[1]] = "C";
	nucl_translation[encoding[2]] = "G";
	nucl_translation[encoding[3]] = "T";
	
	// Lookup translation for Bytes
	for (uint i=0 ; i<256 ; i++) {
		lookup[i] = "";

		for (uint pos=0 ; pos<4 ; pos++) {
			// Get nucleotide
			uint8_t letter = (i >> (2*pos)) & 0b11;
			// Write in the lookup table
			lookup[i] = nucl_translation[letter] + lookup[i];
		}
	}
}


string Stringifyer::translate(const uint8_t * sequence, const size_t nucl_length) const {
	string result = "";
	uint byte_length = nucl_length % 4 == 0 ? nucl_length / 4 : nucl_length / 4 + 1;

	// Prefix can be truncated
	result = lookup[sequence[0]];
	if (nucl_length % 4 != 0) {
		uint to_remove = 4 - (nucl_length % 4);
		result = result.substr(to_remove);
	}

	// Translate Byte per Byte
	for (size_t idx=1 ; idx<byte_length ; idx++) {
		// Translate a Byte
		uint8_t byte = sequence[idx];
		result += lookup[byte];
	}

	return result;
}


string Stringifyer::translate(uint64_t sequence, const size_t nucl_length) const {
	uint8_t bytes[8] = {0, 0, 0, 0, 0, 0, 0, 0};
	for (uint length=0 ; length*4<nucl_length ;	 length++) {
		bytes[7 - length] = sequence & 0xFF;
		sequence >>= 8;
	}

	uint64_t nb_bytes = (nucl_length + 3) / 4;

	return this->translate(bytes + 8 - nb_bytes, nucl_length);
}



Binarizer::Binarizer(const uint8_t encoding[4]) {
	for (uint pos=0 ; pos<4 ; pos++) {
		this->multi_lookup[pos]['A'] = (encoding[0] & 0b11) << (6 - 2 * pos);
		this->multi_lookup[pos]['a'] = (encoding[0] & 0b11) << (6 - 2 * pos);
		this->multi_lookup[pos]['C'] = (encoding[1] & 0b11) << (6 - 2 * pos);
		this->multi_lookup[pos]['c'] = (encoding[1] & 0b11) << (6 - 2 * pos);
		this->multi_lookup[pos]['G'] = (encoding[2] & 0b11) << (6 - 2 * pos);
		this->multi_lookup[pos]['g'] = (encoding[2] & 0b11) << (6 - 2 * pos);
		this->multi_lookup[pos]['T'] = (encoding[3] & 0b11) << (6 - 2 * pos);
		this->multi_lookup[pos]['t'] = (encoding[3] & 0b11) << (6 - 2 * pos);
	}
}


void Binarizer::translate(std::string sequence, uint seq_size, uint8_t * binarized) {
	uint bytes_needed = (seq_size + 3) / 4;

	// First Byte
	binarized[0] = 0;
	for (uint p=(4-(seq_size%4))%4, n=0 ; p<4 ; p++, n++) {
		binarized[0] |= this->multi_lookup[p][sequence[n]];
	}

	// Following bytes
	uint offset = ((seq_size-1) % 4) + 1;
	for (uint b=1 ; b<bytes_needed ; b++) {
		binarized[b] = this->multi_lookup[0][sequence[offset + 4 * (b-1)]];
		binarized[b] |= this->multi_lookup[1][sequence[offset + 4 * (b-1) + 1]];
		binarized[b] |= this->multi_lookup[2][sequence[offset + 4 * (b-1) + 2]];
		binarized[b] |= this->multi_lookup[3][sequence[offset + 4 * (b-1) + 3]];
	}
}

