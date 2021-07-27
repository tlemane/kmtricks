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

namespace km {

struct count_options : km_options
{
  std::string id;
  uint32_t c_ab_min;
  int32_t partition_id;

  std::vector<uint32_t> partition_vec;

  bool clear;
  bool lz4;
  bool kff;
  bool hist;

  std::string format;

  std::string display()
  {
    std::stringstream ss;
    ss << this->global_display();
    RECORD(ss, id);
    RECORD(ss, c_ab_min);
    RECORD(ss, partition_id);
    RECORD(ss, format);
    RECORD(ss, clear);
    RECORD(ss, lz4);
    RECORD(ss, kff);
    RECORD(ss, hist);
    std::string ret = ss.str(); ret.pop_back(); ret.pop_back();
    return ret;
  }
};

using count_options_t = std::shared_ptr<struct count_options>;

};