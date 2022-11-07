#pragma once

#include <vector>
#include <cmath>
#include <bitpacker/bitpacker.hpp>

namespace km {

  std::uint32_t byte_count_pack(std::uint32_t n, std::uint32_t bit)
  {
    return (n * bit + 7) >> 3;
  }

  std::uint32_t to_n_b(std::uint32_t c, std::uint32_t n)
  {
    if (i)
    {
        std::uint64_t r = std::numeric_limits<decltype(i)>::digits - __builtin_clzll(i);

        return r > ((1 << n) - 1) ? ((1 << n) - 1) : r;
    }
    return 0;
  }

  void pack_v(const std::vector<std::uint32_t>& vc, std::vector<uint8_t>& v, int w)
  {
    for (std::size_t i = 0, j = 0; i < vc.size(); ++i, j+=w)
    {
      bitpacker::insert(v, j, w, to_n_b(vc[i], w));
    }
  }
}
