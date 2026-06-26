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
#include <string>
#include <kmtricks/minimizer.hpp>
#include <kmtricks/utils.hpp>

#include <xxhash.h>

namespace km {

class Repartition
{
  inline static const uint32_t s_gatb_magic = 0x12345678;
  Repartition(std::size_t nb_parts, std::size_t nb_minims)
    : m_nb_part(nb_parts), m_nb_minims(nb_minims), m_nb_pass(1), m_has_freq(false)
  {
    m_repart_table.resize(m_nb_minims);
  }

public:
  Repartition(const std::string& path, const std::string& fpath = "")
    : m_path(path), m_fpath(fpath)
  {
    load();
  }

  static Repartition from_xxh(std::size_t nb_partitions, std::size_t minim_size)
  {
    std::size_t nb_minims = std::pow(4, minim_size);
    Repartition repart(nb_partitions, nb_minims);

    for (std::uint32_t m = 0; m < nb_minims; ++m)
    {
      repart.m_repart_table[m] = XXH64(&m, sizeof(m), 0) % nb_partitions;
    }

    return repart;
  }

  void save(const std::string& path) const
  {
    std::ofstream out(path, std::ios::binary | std::ios::out); check_fstream_good(path, out);
    out.write((const char*)&m_nb_part, sizeof(m_nb_part));
    out.write((const char*)&m_nb_minims, sizeof(m_nb_minims));
    out.write((const char*)&m_nb_pass, sizeof(m_nb_pass));
    out.write((const char*)m_repart_table.data(), sizeof(uint16_t)*m_nb_minims);
    out.write((const char*)&m_has_freq, sizeof(m_has_freq));
    out.write((const char*)&s_gatb_magic, sizeof(s_gatb_magic));
  }

  void load()
  {
    std::ifstream in(m_path, std::ios::binary | std::ios::in); check_fstream_good(m_path, in);
    in.read((char*)&m_nb_part, sizeof(m_nb_part));
    in.read((char*)&m_nb_minims, sizeof(m_nb_minims));
    in.read((char*)&m_nb_pass, sizeof(m_nb_pass));
    m_repart_table.resize(m_nb_minims);
    in.read((char*)m_repart_table.data(), sizeof(uint16_t)*m_nb_minims);

    in.read((char*)&m_has_freq, sizeof(m_has_freq));
    in.read((char*)&m_magic, sizeof(m_magic));
    if (m_magic != s_gatb_magic)
      throw IOError("Invalid file format");

    if (m_has_freq && !m_fpath.empty())
    {
      std::ifstream inf(m_fpath, std::ios::binary | std::ios::in);
      m_freq_table.resize(m_nb_minims);
      inf.read(reinterpret_cast<char*>(m_freq_table.data()), sizeof(uint32_t) * m_nb_minims);
      inf.read(reinterpret_cast<char*>(&m_magic), sizeof(m_magic));
      if (m_magic != s_gatb_magic)
        throw IOError("Invalid file format");
    }
  }

  template<size_t MAX_K>
  uint16_t get_partition(const Minimizer<MAX_K>& minim) const
  {
    return m_repart_table[minim.value()];
  }

  uint16_t get_partition(uint32_t value) const
  {
    return m_repart_table[value];
  }

  template<size_t MAX_K>
  uint16_t get_freq_order(const Minimizer<MAX_K>& minim) const
  {
    return m_freq_table[minim.value()];
  }

  uint16_t get_nb_minimizers() const
  {
    return m_nb_minims;
  }

  void write_minimizers(const std::vector<std::string>& paths, size_t size)
  {
    std::vector<std::ofstream> outs;
    for (auto& p: paths)
      outs.push_back(std::ofstream(p, std::ios::out));

    for (size_t i=0; i<m_repart_table.size(); i++)
      outs[m_repart_table[i]] << Mmer(i, size).to_string() << "\n";
  }

  const std::vector<uint16_t>& table() const
  {
    return m_repart_table;
  }

private:
  std::string m_path;
  std::string m_fpath;

  uint16_t m_nb_part;
  uint64_t m_nb_minims;
  uint16_t m_nb_pass;
  bool m_has_freq;
  uint32_t m_magic;
  std::vector<uint16_t> m_repart_table;
  std::vector<uint32_t> m_freq_table;
};

};
