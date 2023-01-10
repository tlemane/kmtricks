#pragma once

#include <bitpacker/bitpacker.hpp>
#include <cmath>
#include <limits>
#include <vector>

namespace km {

/**
 * Compute the number of byte needed to store `n` elements of `bit` bits each.
 * @param n [IN] The number of element
 * @param bit [IN] The size in bit of each element
 */
inline std::uint32_t byte_count_pack(std::uint32_t n, std::uint32_t bit) {
    return (n * bit + 7) >> 3;
}

/**
 * Compute the number of bit needed to represent `c`, caped at `2**n-1`
 * @param c [IN] The element to represent
 * @param n [IN] The maximal size in bit allowed to represent
 */
inline std::uint32_t to_n_b(std::uint32_t c, std::uint32_t n) {
    if (c) {
        // minimum number of bits required to represent c
        std::uint64_t r = std::numeric_limits<decltype(c)>::digits - __builtin_clzll(c);

        // min(r,  2**(n-1))
        return r > n ? ((1 << n) - 1) : c;
    }
    return 0;
}

inline void pack_v(const std::vector<std::uint32_t>& vc, std::vector<uint8_t>& v, int w) {
    for (std::size_t i = 0, j = 0; i < vc.size(); ++i, j += w) {
        bitpacker::insert(v, j, w, to_n_b(vc[i], w));
    }
}
}  // namespace km
