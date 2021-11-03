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

#include <vector>
#include <algorithm>
#include <cstdint>
#include <memory>
#include <iostream>
#include <mutex>

#include <kmtricks/utils.hpp>

namespace km {

enum class KHistType { UNIQUE, TOTAL };

class KHist
{
  template<size_t buf_size>
  friend class HistWriter;

  template<size_t buf_size>
  friend class HistReader;

public:
  KHist() = default;
  KHist(int idx, size_t ksize, size_t lower, size_t upper)
    : m_idx(idx), m_ksize(ksize), m_lower(lower), m_upper(upper)
  {
    m_hist_u.resize(m_upper - m_lower + 1, 0);
    m_hist_n.resize(m_upper - m_lower + 1, 0);
  }

  void inc(uint64_t count)
  {
    m_uniq++;
    m_total += count;
    if (count < m_lower)
    {
      m_oob_lu++;
      m_oob_ln += count;
    }
    else if (count > m_upper)
    {
      m_oob_uu++;
      m_oob_un += count;
    }
    else
    {
      m_hist_u[count - m_lower]++;
      m_hist_n[count - m_lower] += count;
    }
  }

  void set_type(KHistType type)
  {
    m_type = type;
  }

  std::shared_ptr<KHist> clone()
  {
    m_clone = true;
    std::shared_ptr<KHist> hist = std::make_shared<KHist>(m_idx, m_ksize, m_lower, m_upper);
    m_clones.push_back(hist);
    return hist;
  }

  uint64_t unique() const { return m_uniq; }
  uint64_t total() const { return m_total; }
  uint64_t lower() const { return m_lower; }
  uint64_t upper() const { return m_upper; }

  uint64_t oob_lower_unique() const { return m_oob_lu; }
  uint64_t oob_upper_unique() const { return m_oob_uu; }
  uint64_t oob_lower_total() const { return m_oob_ln; }
  uint64_t oob_upper_total() const { return m_oob_un; }

  uint32_t kmer_size() const { return m_ksize; }
  uint32_t idx() const { return m_idx; }

  int64_t get_count(size_t c, KHistType type) const
  {
    if ((c < m_lower) || (c > m_upper))
      return -1;
    if (type == KHistType::UNIQUE)
      return m_hist_u[c];
    return m_hist_n[c];
  }

  const std::vector<uint64_t>& get_vec(KHistType type = KHistType::UNIQUE) const
  {
    if (type == KHistType::UNIQUE)
      return m_hist_u;
    return m_hist_n;
  }

  void merge_clones()
  {
    if (m_clone && !m_merged)
    {
      for (auto& h : m_clones)
      {
        m_uniq += h->m_uniq;
        m_total += h->m_total;
        m_oob_lu += h->m_oob_lu;
        m_oob_uu += h->m_oob_uu;
        m_oob_ln += h->m_oob_ln;
        m_oob_un += h->m_oob_un;
        for (size_t i=0; i<h->m_hist_u.size(); i++)
        {
          m_hist_u[i] += h->m_hist_u[i];
          m_hist_n[i] += h->m_hist_n[i];
        }
      }
      m_merged = true;
      clear_clones();
    }
  }

  void clear_clones()
  {
    m_clones.clear();
    m_clone = false;
    m_merged = false;
  }

  std::string as_string(KHistType type = KHistType::UNIQUE, const std::string sep = "\n") const
  {
    std::stringstream ss;
    uint64_t count = 0;
    auto vec = m_type == KHistType::UNIQUE ? m_hist_u : m_hist_n;
    std::for_each(vec.begin(), vec.end(), [&count, &ss, &sep](uint64_t c){
      ss << std::to_string(count) << " " << std::to_string(c) << sep;
      count++;
    });
    return ss.str();
  }

  void print(std::ostream& output_stream, KHistType type = KHistType::UNIQUE, const std::string sep = "\n") const
  {
    output_stream << as_string(type, sep);
  }

  auto begin()
  {
    if (m_type == KHistType::UNIQUE)
      return m_hist_u.begin();
    return m_hist_n.begin();
  }

  auto end()
  {
    if (m_type == KHistType::UNIQUE)
      return m_hist_u.end();
    return m_hist_n.end();
  }

  auto cbegin() const
  {
    if (m_type == KHistType::UNIQUE)
      return m_hist_u.cbegin();
    return m_hist_n.cbegin();
  }

  auto cend()
  {
    if (m_type == KHistType::UNIQUE)
      return m_hist_u.cend();
    return m_hist_n.cend();
  }

private:
  int32_t m_idx {0};
  uint32_t m_ksize {0};
  uint64_t m_lower {0};
  uint64_t m_upper {0};
  uint64_t m_uniq {0};
  uint64_t m_total {0};
  uint64_t m_oob_lu {0};
  uint64_t m_oob_uu {0};
  uint64_t m_oob_ln {0};
  uint64_t m_oob_un {0};

  std::vector<uint64_t> m_hist_u;
  std::vector<uint64_t> m_hist_n;
  std::vector<std::shared_ptr<KHist>> m_clones;
  KHistType m_type {KHistType::UNIQUE};
  bool m_clone {false};
  bool m_merged {false};
};

using hist_t = std::shared_ptr<KHist>;

inline hist_t get_hist_clone(hist_t hist)
{
  if (hist)
    return hist->clone();
  return nullptr;
}

inline std::vector<uint32_t> compute_merge_thresholds(std::vector<hist_t>& histograms,
                                                      double p,
                                                      const std::string& path)
{
  std::vector<uint32_t> thresholds(histograms.size());
  for (size_t h=0; h<histograms.size(); h++)
  {
    uint32_t sum = 0;
    uint32_t n = histograms[h]->unique() * p;
    auto v = histograms[h]->get_vec(KHistType::UNIQUE);
    for (size_t i=0; i<v.size(); i++)
    {
      if (sum > n)
      {
        thresholds.push_back(i);
        break;
      }
      sum += v[i];
    }
  }
  std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
  for (auto& t: thresholds)
  {
    out << std::to_string(t) << "\n";
  }
  return thresholds;
}

};
