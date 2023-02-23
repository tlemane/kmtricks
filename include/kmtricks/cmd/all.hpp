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
#include <kmtricks/exceptions.hpp>
#include <kmtricks/utils.hpp>
#include <kmtricks/io/fof.hpp>

namespace km {

struct all_options : km_options
{
  std::string fof;
  uint32_t kmer_size {0};
  uint32_t c_ab_min {0};
  uint32_t m_ab_min {0};
  uint32_t r_min {0};
  std::string m_ab_min_path;
  double m_ab_min_f {0.0};
  bool m_ab_float = {false};
  uint32_t save_if {0};

  uint32_t minim_type {0};
  uint32_t minim_size {0};
  uint32_t repart_type {0};
  uint32_t nb_parts {0};

  uint64_t bloom_size {0};

  bool keep_tmp {false};
  bool lz4 {false};
  bool kff {false};
  bool skip_merge {false};
  bool hist {false};

  uint32_t bwidth {0};

  uint32_t max_memory {8000};
  double restrict_to;
  std::vector<uint32_t> restrict_to_list;
  std::vector<uint32_t> m_ab_min_vec;

  double focus {1.0};

  std::string from;

  MODE mode;
  FORMAT format;
  OUT_FORMAT out_format;
  COUNT_FORMAT count_format;
  COMMAND until;

#ifdef WITH_PLUGIN
  std::string plugin;
  std::string plugin_config;
  bool use_plugin {false};
#endif

  std::string display()
  {
    std::stringstream ss;
    ss << this->global_display();
    RECORD(ss, fof);
    RECORD(ss, kmer_size);
    RECORD(ss, c_ab_min);
    RECORD(ss, m_ab_min);
    RECORD(ss, r_min);
    RECORD(ss, m_ab_min_path);
    RECORD(ss, m_ab_min_f);
    RECORD(ss, m_ab_float);
    RECORD(ss, save_if);
    RECORD(ss, minim_size);
    RECORD(ss, minim_type);
    RECORD(ss, repart_type);
    RECORD(ss, nb_parts);
    RECORD(ss, bloom_size);
    RECORD(ss, keep_tmp);
    RECORD(ss, lz4);
    RECORD(ss, kff);
    RECORD(ss, skip_merge);
    RECORD(ss, hist);
    RECORD(ss, focus);
    RECORD(ss, restrict_to);
    RECORD(ss, bwidth);
#ifdef WITH_PLUGIN
    RECORD(ss, use_plugin);
    RECORD(ss, plugin);
    RECORD(ss, plugin_config);
#endif
    ss << "mode=" << mode_to_str(mode) << ", ";
    ss << "format=" << format_to_str2(format) << ", ";
    ss << "bf_format=" << format_to_str(out_format) << ", ";
    ss << "count_format=" << cformat_to_str(count_format) << ", ";
    ss << "until=" << cmd_to_str(until);
    return ss.str();
  }

  void sanity_check()
  {
    if ((kff) && (until != COMMAND::COUNT))
    {
      throw PipelineError("--kff-output/--kff-sk-output available only with --until count");
    }
    if ((kff) && (count_format == COUNT_FORMAT::HASH))
    {
      throw PipelineError("--kff-output/--kff-sk-output available only in k-mer mode.");
    }
    if (skip_merge)
    {
      if ((mode != MODE::BFT) || (count_format != COUNT_FORMAT::HASH))
      {
        throw PipelineError("--skip-merge available only with --mode hash:bft:bin");
      }
    }
    if ((mode == MODE::BFT || mode == MODE::BF))
    {
      if ((restrict_to != 1.0) || !restrict_to_list.empty())
      {
        throw PipelineError("--mode bf|bft requires all partitions.");
      }
    }
    Fof fof_file(fof);

    if (m_ab_float)
    {
      hist = true;
    }
    else if (!m_ab_min_path.empty())
    {
      std::ifstream in(m_ab_min_path, std::ios::in); check_fstream_good(m_ab_min_path, in);
      for (std::string line; std::getline(in, line);)
      {
        m_ab_min_vec.push_back(std::stol(line));
      }
      if (fof_file.size() != m_ab_min_vec.size())
      {
        throw PipelineError(fmt::format(
          "The number of thresholds in {} is different from the number of samples.", m_ab_min_path));
      }
    }
    else
    {
      for (size_t i=0; i<fof_file.size(); i++)
      {
        m_ab_min_vec.push_back(m_ab_min);
      }
    }
  }

  void dump(const std::string& path)
  {
    std::ofstream out_opt(path, std::ios::out); check_fstream_good(path, out_opt);
    out_opt << display();
  }
};

using all_options_t = std::shared_ptr<struct all_options>;

};
