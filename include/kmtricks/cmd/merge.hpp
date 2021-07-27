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
#include <memory>
#include <thread>

#include <spdlog/spdlog.h>


#include <kmtricks/cli/cli_common.hpp>
#include <kmtricks/cmd/cmd_common.hpp>
#include <kmtricks/kmdir.hpp>

namespace km {

struct merge_options : km_options
{
  uint32_t m_ab_min;
  std::string m_ab_min_path;
  double m_ab_min_f;
  bool m_ab_float;
  uint32_t r_min;
  int32_t partition_id;
  uint32_t save_if;
  std::vector<uint32_t> m_ab_min_vec;

  bool clear;
  bool lz4;

  MODE mode;
  FORMAT format;
  COUNT_FORMAT count_format;

  std::string display()
  {
    std::stringstream ss;
    ss << this->global_display();
    RECORD(ss, m_ab_min);
    RECORD(ss, m_ab_min_path);
    RECORD(ss, r_min);
    RECORD(ss, partition_id);
    RECORD(ss, save_if);
    RECORD(ss, clear);
    RECORD(ss, lz4);
    std::string ret = ss.str(); ret.pop_back(); ret.pop_back();
    return ret;
  }

  void init_vector()
  {
    if (!m_ab_min_path.empty())
    {
      std::ifstream in(m_ab_min_path, std::ios::in); check_fstream_good(m_ab_min_path, in);
      for (std::string line; std::getline(in, line);)
      {
        m_ab_min_vec.push_back(std::stol(line));
      }
      if (KmDir::get().m_fof.size() != m_ab_min_vec.size())
      {
        throw PipelineError(fmt::format(
          "The number of thresholds in {} is different from the number of samples.", m_ab_min_path));
      }
    }
    else
    {
      for (size_t i=0; i<KmDir::get().m_fof.size(); i++)
      {
        m_ab_min_vec.push_back(m_ab_min);
      }
    }
  }
};

using merge_options_t = std::shared_ptr<struct merge_options>;

};