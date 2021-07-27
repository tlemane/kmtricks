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

class KHist
{
public:
  KHist() {}
  KHist(int idx, size_t ksize, size_t lower, size_t upper)
    : idx(idx), ksize(ksize), lower(lower), upper(upper)
  {
    hist_u.resize(upper-lower+1, 0);
    hist_n.resize(upper-lower+1, 0);
  }

  void inc(uint64_t count)
  {
    std::unique_lock<std::mutex> lock(m_mutex);
    uniq++;
    total += count;
    if (count < lower)
    {
      oob_lu++;
      oob_ln+=count;
    }
    else if (count > upper)
    {
      oob_uu++;
      oob_un+=count;
    }
    else
    {
      hist_u[count-lower]++;
      hist_n[count-lower]+=count;
    }
  }

  void print_histu()
  {
    _print(hist_u);
  }

  void print_histn()
  {
    _print(hist_n);
  }

private:
  void _print(const std::vector<uint64_t>& v)
  {
    uint64_t current = lower;
    for_each(v.begin(), v.end(), [&current](uint64_t c){
      std::cerr << std::to_string(current) << " " << std::to_string(c) << "\n";
      current++;
    });
    std::cerr << std::flush;
  }

public:
  int32_t  idx {0};
  uint32_t ksize {0};
  uint64_t lower {0};
  uint64_t upper {0};
  uint64_t uniq {0};
  uint64_t total {0};
  uint64_t oob_lu {0};
  uint64_t oob_uu {0};
  uint64_t oob_ln {0};
  uint64_t oob_un {0};
  std::vector<uint64_t> hist_u;
  std::vector<uint64_t> hist_n;
  std::mutex m_mutex;
};

using hist_t = std::shared_ptr<KHist>;

inline std::vector<uint32_t> compute_merge_thresholds(std::vector<hist_t>& histograms,
                                                      double p,
                                                      const std::string& path)
{
  std::vector<uint32_t> thresholds(histograms.size());
  for (size_t h=0; h<histograms.size(); h++)
  {
    uint32_t sum = 0;
    uint32_t n = histograms[h]->uniq * p;
    for (size_t i=0; i<histograms[h]->hist_u.size(); i++)
    {
      if (sum > n)
      {
        thresholds.push_back(i);
        break;
      }
      sum += histograms[h]->hist_u[i];
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