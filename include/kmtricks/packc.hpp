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
    if (c == 1)
      return 1;

    if (c == 2)
      return 2;

    auto c_ = std::ceil(std::log2(c));
    auto p = static_cast<std::uint32_t>(std::pow(2, n)) - 1;

    return c_ > p ? p : c_;
  }

  void pack_v(const std::vector<std::uint32_t>& vc, std::vector<uint8_t>& v, int w)
  {
    for (std::size_t i = 0, j = 0; i < vc.size(); ++i, j+=w)
    {
      spdlog::debug("{} -> {}", vc[i], to_n_b(vc[i], w));
      bitpacker::insert(v, j, w, to_n_b(vc[i], w));
    }
  }
}
