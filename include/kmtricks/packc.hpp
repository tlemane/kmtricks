#pragma once

#include <bitpacker/bitpacker.hpp>
#include <cmath>
#include <limits>
#include <vector>

namespace km
{

/**
 * Compute the number of byte needed to store `n` elements of `bit` bits each.
 * @param n [IN] The number of element
 * @param bit [IN] The size in bit of each element
 */
inline std::uint32_t byte_count_pack(std::uint32_t n, std::uint32_t bit)
{
  return (n * bit + 7) >> 3;
}

/**
 * Compute ceil(log2(`c` + 1)), caped at `2**n-1`
 * @param c [IN] The argument of log2
 * @param max_width [IN] maximal width (in bits) of the image of log2
 */
inline std::uint32_t to_n_b(std::uint32_t c, std::uint32_t max_width)
{
  if (c)
  {
    std::uint64_t r = std::numeric_limits<decltype(c)>::digits - __builtin_clz(c);
    return r > ((1 << max_width) - 1) ? ((1 << max_width) - 1) : r;
  }
  return 0;
}

template<typename C>
inline void pack_v(const std::vector<C>& vc, std::vector<uint8_t>& v, int w)
{
  for (std::size_t i = 0, j = 0; i < vc.size(); ++i, j += w)
  {
    bitpacker::insert(v, j, w, to_n_b(vc[i], w));
  }
}
}  // namespace km
