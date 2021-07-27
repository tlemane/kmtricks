// hash function from HowDeSBT, B. Harris, P. Medvedev
// (https://github.com/medvedevgroup/HowDeSBT) 

// sabuhash.h-- implementation of a variant of nthash [1], replacing rol() with
//              an xor-shift generator.

// The hash function is implemented in both canonical and non-canonical
// versions, and rolling and non-rolling forms. The rolling form is expected to
// be faster, especially as k (the word size) increases. The caller can choose
// whether or not N should be considered a valid character.
//
// The basic hash function is seed-based, the intent being that different seeds
// can be used to create as many hash functions as needed. However, if more
// than two hash functions are needed, it is more efficient to generate
// additional hash values from two basis functions, using fill_hash_values().
//
// The hash function reports zero to indicate an invalid result (e.g. when the
// kmer includes an invalid character). On the rare occasion when the hash
// function value mathematically should be zero, the value one is reported
// instead.
//
// Implementation details:
//
// The core of the hash is, like ntHash [1], the 64-bit xor-sum of position-
// and nucleotide-specific values. In contrast to ntHash, an xor-shift type
// function forward() is used to account for a nucleotide's position (replacing
// ntHash's use of rotate left). This eliminates the possibility of two like
// nucleotides 64 bp apart cancelling each other.
//
// An avalanche step, borrowed from murmurhash3 [4a], is applied to the core
// result so that the resulting hash function has better randomness properties.
//
// Mathematically, we have
//   h = a(xor {for i=0..k-1} of f(z[s[k-1-i]],i))
// where
//   s is the kmer
//   k is the length of the kmer
//   s[i] is the nucleotide in position i of the kmer
//   z[c] is the kernel value for nucleotide c (one of kernelA, kernelC, ...)
//   f(u,i) is the application of forward() to value u, i times (i may be zero)
//   a(u) is the application of avalanche() to value u
//   h is the hash value of s
//
// In other words,
//   u = f(z[s[0]],k-1)) xor ... xor f(z[s[k-2]],1)) xor z[s[k-1]]
//   h = a(u)
//
// The xor-shift function forward(u) is implemented in two steps
//   u = u ^ (u << L)
//   u = u ^ (u >> R)
// where L and R have been chosen such that repeated application of forward()
// doesn't have (many) short-period cycles, and is easily invertible. Ease of
// invertibility is accomplished if 16<L<32<R<64.
//
// The five kernel values (for A, C, G, T and N) were randomly chosen and
// satisfy certain properties. In particular, they are all separated by at
// least 100 billion steps of the forward() function, have period at least 100
// billion, and none of the first 500 thousand steps is heavily skewed to 0s or
// 1s (all have at least 16 0s and 16 1s). Since the hash function core is a
// mix of some subset of the first k values of forward() starting with any of
// the four (or five) kernels, the goal was that these 4k (or 5k) values should
// be distinct and all contribute both 0s and 1s to the xor result.
//
// The canonical version computes the non-canonical core hash of both the
// forward and reverse-complement kmer, and returns their sum (prior to the
// avalanche step). Note that many canonical hash implementations use
// min(forward,reverse), but that biases the result toward lower values. (That
// bias, though, would probably be obscured here by the avalanche step). Also
// note that at least one canonical hash implementation uses
// xor(forward,reverse), but that causes collisions for all palindromic kmers.
// (And that fact can't be corrected by the avalanche step).
//
// When the hash function is seeded with 63-bit seed s, the kernel values are
// multiplied by (2s+1). Multiplying by an odd value guarantees than these
// new kernels will remain non-zero. Hypothetically the resulting function will
// be just as random as (and uncorrelated with) one with a different seed, but
// at present this has not been well tested. The 100 billion step separation
// and period properties have been tested for seeds 0, 1, 2, 3, and 4. The 500
// thousand step non-skew property has only been tested for seed 0.
//
// fill_hash_values() uses values from two different hash functions as seeds
// for a PRNG to generate a list of hash values. This idea and the choice of
// PRNG are suggested by the implementation of BBHash [2]. The PRNG is
// xorshift128+ due to Vigna, described in [3a], implemented at [3b].
//
// References:
//
//   [1]  Mohamadi, Hamid, et al. "ntHash: recursive nucleotide hashing."
//        Bioinformatics 32.22 (2016): 3492-3494.
//   [2]  Limasset, Antoine, et al. "Fast and scalable minimal perfect hashing
//        for massive key sets." arXiv preprint arXiv:1702.03154 (2017).
//   [3a] Vigna, Sebastiano. "Further scramblings of Marsaglia’s xorshift
//        generators." Journal of Computational and Applied Mathematics 315
//        (2017): 175-181.
//   [3b] https://github.com/jj1bdx/xorshiftplus
//        https://github.com/jj1bdx/xorshiftplus/blob/master/c/xorshift128plus-speed.c
//   [4a] Austin Appleby. Smhasher and murmurhash3 webpage, 2016.
//        https://github.com/ aappleby/smhasher/wiki
//   [4b] https://github.com/aappleby/smhasher/wiki/MurmurHash3

#ifndef sabuhash_H
#define sabuhash_H

#include <cassert>
#include <cstdint>
#include <string>

#define sabuhash_xorShiftL 17 // 16 < xorShiftL < 32
#define sabuhash_xorShiftR 47 // 32 < xorShiftR < 64
#define sabuhash_kernelA 0xC02069C411718AC9
#define sabuhash_kernelC 0x29A37009B8869707
#define sabuhash_kernelG 0x1A14AE384351C3F6
#define sabuhash_kernelT 0x900E988F0E40231E
#define sabuhash_kernelN 0x8C41F1213A951881

#define sabuhash_Abits 0
#define sabuhash_Cbits 1
#define sabuhash_Gbits 2
#define sabuhash_Tbits 3

class SabuHash
{
public:
	static const int xorShiftL = sabuhash_xorShiftL;
	static const int xorShiftR = sabuhash_xorShiftR;
	static const std::uint64_t kernelA = sabuhash_kernelA;
	static const std::uint64_t kernelC = sabuhash_kernelC;
	static const std::uint64_t kernelG = sabuhash_kernelG;
	static const std::uint64_t kernelT = sabuhash_kernelT;
	static const std::uint64_t kernelN = sabuhash_kernelN;

	unsigned int k;
	std::uint64_t seed;
	bool allowN;
	bool validNt; // true => latest nucleotide was valid
	unsigned int charsAccumulated;
	std::uint64_t kernelTable[256], forwardByK[256], forwardByK2[4];
	std::uint64_t hForward;

public:
	SabuHash(const unsigned int _k,
					 const std::uint64_t _seed = 0,
					 const bool _allowN = false)
			: k(_k),
				seed(_seed),
				allowN(_allowN),
				validNt(false),
				charsAccumulated(0),
				hForward(0)
	{
		assert(k > 0);

		std::uint64_t initialA, initialC, initialG, initialT, initialN;
		std::uint64_t forwardA, forwardC, forwardG, forwardT, forwardN;

		forwardA = initialA = kernelA * (2 * _seed + 1);
		forwardC = initialC = kernelC * (2 * _seed + 1);
		forwardG = initialG = kernelG * (2 * _seed + 1);
		forwardT = initialT = kernelT * (2 * _seed + 1);
		forwardN = initialN = kernelN * (2 * _seed + 1);
		for (unsigned int i = 0; i < k; i++)
		{
			forwardA = forward(forwardA);
			forwardC = forward(forwardC);
			forwardG = forward(forwardG);
			forwardT = forward(forwardT);
			forwardN = forward(forwardN);
		}

		for (unsigned int ch = 0; ch < 256; ch++)
			kernelTable[ch] = forwardByK[ch] = 0;

		kernelTable[(unsigned char)'A'] = kernelTable[(unsigned char)'a'] = initialA;
		kernelTable[(unsigned char)'C'] = kernelTable[(unsigned char)'c'] = initialC;
		kernelTable[(unsigned char)'G'] = kernelTable[(unsigned char)'g'] = initialG;
		kernelTable[(unsigned char)'T'] = kernelTable[(unsigned char)'t'] = initialT;
		forwardByK[(unsigned char)'A'] = forwardByK[(unsigned char)'a'] = forwardA;
		forwardByK[(unsigned char)'C'] = forwardByK[(unsigned char)'c'] = forwardC;
		forwardByK[(unsigned char)'G'] = forwardByK[(unsigned char)'g'] = forwardG;
		forwardByK[(unsigned char)'T'] = forwardByK[(unsigned char)'t'] = forwardT;

		if (allowN)
		{
			kernelTable[(unsigned char)'N'] = kernelTable[(unsigned char)'n'] = initialN;
			forwardByK[(unsigned char)'N'] = forwardByK[(unsigned char)'n'] = forwardN;
		}

		forwardByK2[sabuhash_Abits] = forwardA;
		forwardByK2[sabuhash_Cbits] = forwardC;
		forwardByK2[sabuhash_Gbits] = forwardG;
		forwardByK2[sabuhash_Tbits] = forwardT;
	}

	~SabuHash() {}

	inline std::uint64_t hash(const char *s) // only first k characters are used
	{
		// note that we return 0 iff we don't have k valid characters

		std::uint64_t hF = 0;
		for (unsigned int i = 0; i < k; i++)
		{
			unsigned char chIn = s[i];
			if (kernelTable[chIn] == 0)
				return 0;
			hF = forward(hF) ^ kernelTable[chIn];
		}

		hForward = hF;

		hF = avalanche(hF);
		if (hF == 0)
			return 1;
		else
			return hF;
	}

	inline std::uint64_t hash(const std::string &s) // only first k characters are used
	{
		// note that we return 0 iff we don't have k valid characters

		if (s.length() < k)
			return 0;

		std::uint64_t hF = 0;
		for (unsigned int i = 0; i < k; i++)
		{
			unsigned char chIn = s[i];
			if (kernelTable[chIn] == 0)
				return 0;
			hF = forward(hF) ^ kernelTable[chIn];
		}

		hForward = hF;

		hF = avalanche(hF);
		if (hF == 0)
			return 1;
		else
			return hF;
	}

	inline std::uint64_t hash(const std::uint64_t *_data) // kth nuc in lsbits, k-31st nuc in msbits
	{																											// k-32nd nuc in lsbits of next word, etc.
		std::uint64_t *data = (std::uint64_t *)_data;

		// note that we never return 0

		std::uint64_t hF = 0;
		std::uint64_t d = 0;
		for (unsigned int i = 0; i < k; i++)
		{
			if ((i % 32) == 0)
				d = *(data++);
			unsigned char twoBits = d & 3;
			d >>= 2;
			hF = backward(hF ^ forwardByK2[twoBits]);
		}

		hForward = hF;

		hF = avalanche(hF);
		if (hF == 0)
			return 1;
		else
			return hF;
	}

	inline std::uint64_t rolling_hash(const unsigned char chIn,
																		const unsigned char chOut = 0) // this was chIn k characters earlier
	{
		// note that we return 0 iff we don't yet have k characters

		if (kernelTable[chIn] == 0)
		{
			// chIn not in {A,C,G,T}
			validNt = false;
			charsAccumulated = 0;
			hForward = 0;
			return 0;
		}

		validNt = true;
		if (charsAccumulated < k)
		{
			hForward = forward(hForward) ^ kernelTable[chIn];
			charsAccumulated++;
		}
		else
		{
			// we assume, without checking, that chOut is in {A,C,G,T}
			hForward = forward(hForward) ^ kernelTable[chIn] ^ forwardByK[chOut];
		}

		if (charsAccumulated < k)
			return 0;

		std::uint64_t h = avalanche(hForward);
		if (h == 0)
			return 1;
		else
			return h;
	}

	inline void reset_rolling_hash()
	{
		rolling_hash(0);
	}

	static inline std::uint64_t forward(const std::uint64_t v)
	{
		std::uint64_t u = v;
		u ^= u << xorShiftL;
		u ^= u >> xorShiftR;
		return u;
	}

	static inline std::uint64_t backward(const std::uint64_t v)
	{
		// assumes 16 < xorShiftL < 32 < xorShiftR < 64
		std::uint64_t u = v;
		u ^= u >> xorShiftR;
		u ^= u << xorShiftL;
		u ^= u << (2 * xorShiftL);
		return u;
	}

	static inline std::uint64_t avalanche(const std::uint64_t v)
	{
		// 64-bit avalanche function, from murmurhash3; see [4a] and [4b]
		std::uint64_t u = v;
		u ^= u >> 33;
		u *= 0xFF51AFD7ED558CCD;
		u ^= u >> 33;
		u *= 0xC4CEB9FE1A85EC53;
		u ^= u >> 33;
		return u;
	}

	static inline std::uint64_t ehcnalava(const std::uint64_t v)
	{
		// inverse of murmurhash3's 64-bit avalanche function
		// nota bene: 0x9CB4B2F8129337DB*0xC4CEB9FE1A85EC53 == 1 mod 2^64
		//            0x4F74430C22A54005*0xFF51AFD7ED558CCD == 1 mod 2^64
		//            u ^= u >> 33 is its own inverse
		//            "ehcnalava" is "avalanche" spelt backwards
		std::uint64_t u = v;
		u ^= u >> 33;
		u *= 0x9CB4B2F8129337DB;
		u ^= u >> 33;
		u *= 0x4F74430C22A54005;
		u ^= u >> 33;
		return u;
	}

	static inline void fill_hash_values(std::uint64_t hashValues[],
																			const int numHashes,
																			const std::uint64_t h1,
																			const std::uint64_t h2)
	{
		// use two hash values to seed a PRNG, then generate N-2 additional
		// hash values; idea suggested by the implementation of BBHash [2]; the
		// PRNG is xorshift128+, see [3a] and [3b]

		std::uint64_t s0, s1, s2;
		hashValues[0] = s0 = h1;
		hashValues[1] = s1 = h2;
		for (int h = 2; h < numHashes; h++)
		{ // see [3a] and [3b]
			s1 ^= s1 << 23;
			hashValues[h] = s2 = s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26);
			s0 = s1;
			s1 = s2;
		}
	}
};

#undef sabuhash_xorShiftL
#undef sabuhash_xorShiftR
#undef sabuhash_kernelA
#undef sabuhash_kernelC
#undef sabuhash_kernelG
#undef sabuhash_kernelT

#endif // sabuhash_H
