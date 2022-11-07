#pragma once

#include <vector>
#include <cmath>
#include <limits>
#include <bitpacker/bitpacker.hpp>

namespace km {

  inline std::uint32_t byte_count_pack(std::uint32_t n, std::uint32_t bit)
  {
    return (n * bit + 7) >> 3;
  }

  inline std::uint32_t to_n_b(std::uint32_t c, std::uint32_t n)
  {
    if (c)
    {
        std::uint64_t r = std::numeric_limits<decltype(c)>::digits - __builtin_clzll(c);

        return r > ((1 << n) - 1) ? ((1 << n) - 1) : r;
    }
    return 0;
  }

  inline void pack_v(const std::vector<std::uint32_t>& vc, std::vector<uint8_t>& v, int w)
  {
    for (std::size_t i = 0, j = 0; i < vc.size(); ++i, j+=w)
    {
      bitpacker::insert(v, j, w, to_n_b(vc[i], w));
    }
  }
}
