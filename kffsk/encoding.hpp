#include <cstdint>
#include <iostream>
#include <string>
#include <map>


#ifndef ENCODING_H
#define ENCODING_H

/** 
	* Object used to translate compacted sequences from an encoding to another.
	* A translation Byte size lookup table is computed during the object creation
	* and used to fast translate sequences Byte per Byte.
	**/
class Translator {
private:
	uint8_t lookup[256];

public:
	/**
		* Constructor that construct a 256 Bytes lookup table for fast translation
		*
		* @param source The source encoding. A 4 cell array where each of the numbers
		* from 0 to 3 must be present.
		* @param destination The destination encoding. A 4 cell array where each of
		* the numbers from 0 to 3 must be present.
		**/
	Translator(uint8_t source[4], uint8_t destination[4]);
	/**
	  * Inplace translate the sequence from the source encoding to the destination.
	  *
	  * @param sequence 2-bit compacted sequence that will be translated regarding
	  * the encodings.
	  * @param byte_length Length in Bytes of the sequence.
	  **/
	void translate(uint8_t * sequence, size_t byte_length);
};

class RevComp {
public:
	uint8_t reverse[4];
	uint8_t translations[256];

	RevComp(const uint8_t encoding[4]);

	void rev_comp(uint8_t * seq, const uint64_t seq_size) const;
	void rev_data(uint8_t * data, const uint64_t data_size, const uint64_t nb_kmers) const;

	static uint rev_position(const uint fwd_pos, const uint seq_size) {
		return seq_size - fwd_pos - 1;
	}
};

/** 
	* Object used to convert compacted sequences from 2 bits per nucleotide to string.
	* A translation Byte size lookup table is computed during the object creation
	* and used to fast transform sequences Byte per Byte.
	**/
class Stringifyer {
private:
	std::string lookup[256];

public:
	/**
		* Construct a 256 Bytes lookup table for fast conversion
		*
		* @param encoding The sequence encoding
		**/
	Stringifyer(uint8_t encoding[4]);
	/**
	  * read a 2-bits/nucl sequence and return a coresponding string
	  *
	  * @param sequence 2-bit compacted sequence that will be converted regarding
	  * the encodings.
	  * @param nucl_length Length in nucleotides of the sequence.
	  **/
	std::string translate(const uint8_t * sequence, const size_t nucl_length) const;
	std::string translate(uint64_t sequence, const size_t nucl_length) const;
};


class Binarizer {
private:
	std::map<std::string, uint8_t> lookup;
	std::map<char, uint8_t> multi_lookup[4];

public:
	/**
		* Construct a lookup table to translate strings to uint8_t values.
		* This lookup table is used to translate longer than 4 strings.
		*
		* @param encoding The 2-bit nucleotide encoding
		*/
	Binarizer(const uint8_t encoding[4]);
	/**
		* Translate the content of the sequence into a binarized version.
		*
		* @param sequence string sequence to translate.
		* @param binarized array where the binary version is stored.
		* The space must be allocated outside of the function.
		*/
	void translate(std::string sequence, uint seq_size, uint8_t * binarized);
};

#endif
