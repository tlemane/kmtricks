/*****************************************************************************
 *   kmtricks
 *   Authors: T. Lemane
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#pragma once
#include <fstream>
#include <cmath>

#include <kmtricks/utils.hpp>

namespace km {

class HashWindow
{
public:
  HashWindow() {}
  HashWindow(uint64_t bloom_size, uint64_t nb_partitions, uint32_t minimizer_size)
    : m_nb_partitions(nb_partitions), m_minim_size(minimizer_size)
  {
    m_window_size_bits = ROUND_UP(
      static_cast<uint64_t>(
        std::ceil(static_cast<double>(bloom_size)/static_cast<double>(m_nb_partitions))), 64);
    m_window_size_bytes = NBYTES(m_window_size_bits);
    m_bloom_size = m_window_size_bits * nb_partitions;
    assert(m_window_size_bits % 8 == 0);
  }

  HashWindow(const std::string& path)
  {
    std::ifstream in(path, std::ios::binary | std::ios::in); check_fstream_good(path, in);
    in.read(reinterpret_cast<char*>(&m_bloom_size), sizeof(m_bloom_size));
    in.read(reinterpret_cast<char*>(&m_nb_partitions), sizeof(m_nb_partitions));
    in.read(reinterpret_cast<char*>(&m_window_size_bits), sizeof(m_window_size_bits));
    in.read(reinterpret_cast<char*>(&m_window_size_bytes), sizeof(m_window_size_bytes));
    in.read(reinterpret_cast<char*>(&m_minim_size), sizeof(m_minim_size));
  }

  void serialize(const std::string& path)
  {
    std::ofstream out(path, std::ios::binary | std::ios::out); check_fstream_good(path, out);
    out.write(reinterpret_cast<char*>(&m_bloom_size), sizeof(m_bloom_size));
    out.write(reinterpret_cast<char*>(&m_nb_partitions), sizeof(m_nb_partitions));
    out.write(reinterpret_cast<char*>(&m_window_size_bits), sizeof(m_window_size_bits));
    out.write(reinterpret_cast<char*>(&m_window_size_bytes), sizeof(m_window_size_bytes));
    out.write(reinterpret_cast<char*>(&m_minim_size), sizeof(m_minim_size));
  }

  uint64_t get_window_size_bytes() const
  {
    return m_window_size_bytes;
  }

  uint64_t get_window_size_bits() const
  {
    return m_window_size_bits;
  }

  uint64_t get_lower(uint32_t partition_id) const
  {
    return partition_id * m_window_size_bits;
  }

  uint64_t get_upper(uint32_t partition_id) const
  {
    return ((partition_id + 1) * m_window_size_bits) - 1;
  }

  uint64_t bloom_size() const
  {
    return m_bloom_size;
  }

  uint32_t minim_size() const
  {
    return m_minim_size;
  }

private:
  uint64_t m_bloom_size {0};
  uint64_t m_nb_partitions {0};
  uint64_t m_window_size_bits {0};
  uint64_t m_window_size_bytes {0};
  uint32_t m_minim_size {0};
};

};