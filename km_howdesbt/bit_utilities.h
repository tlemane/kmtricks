#ifndef bit_utilities_H
#define bit_utilities_H

//----------
//
// prototypes for functions in this module--
//
//----------

bool          bitwise_is_all_zeros (const void* bits, const std::uint64_t numBits);
bool          bitwise_is_all_ones  (const void* bits, const std::uint64_t numBits);
void          bitwise_copy         (const void* bits, void* dstBits, const std::uint64_t numBits);
void          bitwise_and          (const void* bits1, const void* bits2, void* dstBits, const std::uint64_t numBits);
void          bitwise_and          (void* dstBits, const void* bits2, const std::uint64_t numBits);
void          bitwise_mask         (const void* bits1, const void* maskBits, void* dstBits, const std::uint64_t numBits);
void          bitwise_mask         (void* dstBits, const void* maskBits, const std::uint64_t numBits);
void          bitwise_or           (const void* bits1, const void* bits2, void* dstBits, const std::uint64_t numBits);
void          bitwise_or           (void* dstBits, const void* bits2, const std::uint64_t numBits);
void          bitwise_or_not       (const void* bits1, const void* bits2, void* dstBits, const std::uint64_t numBits);
void          bitwise_or_not       (void* dstBits, const void* bits2, const std::uint64_t numBits);
void          bitwise_xor          (const void* bits1, const void* bits2, void* dstBits, const std::uint64_t numBits);
void          bitwise_xor          (void* dstBits, const void* bits2, const std::uint64_t numBits);
void          bitwise_xnor         (const void* bits1, const void* bits2, void* dstBits, const std::uint64_t numBits);
void          bitwise_xnor         (void* dstBits, const void* bits2, const std::uint64_t numBits);
void          bitwise_complement   (const void* bits, void* dstBits, const std::uint64_t numBits);
void          bitwise_complement   (void* dstBits, const std::uint64_t numBits);
void          bitwise_fill         (void* dstBits, const int bitVal, const std::uint64_t numBits);
std::uint64_t bitwise_squeeze      (const void* bits, const void* specBits, const std::uint64_t numBits, void* dstBits,
                                    const std::uint64_t numDstBits=((std::uint64_t)-1));
std::uint64_t bitwise_unsqueeze    (const void* bits, const std::uint64_t numBits,
                                    const void* specBits, const std::uint64_t numSpecBits,
                                    void* dstBits, const std::uint64_t numDstBits=((std::uint64_t)-1));
std::uint64_t bitwise_count        (const void* bits, const std::uint64_t numBits);
std::uint64_t bitwise_and_count    (const void* bits1, const void* bits2, const std::uint64_t numBits);
std::uint64_t bitwise_mask_count   (const void* bits1, const void* bits2, const std::uint64_t numBits);
std::uint64_t bitwise_or_count     (const void* bits1, const void* bits2, const std::uint64_t numBits);
std::uint64_t bitwise_or_not_count (const void* bits1, const void* bits2, const std::uint64_t numBits);
std::uint64_t bitwise_xor_count    (const void* bits1, const void* bits2, const std::uint64_t numBits);

#define hamming_distance(bits1,bits2,numBits) (bitwise_xor_count((bits1),(bits2),(numBits)))

#endif // bit_utilities_H
