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

struct agg_options : km_options
{
  std::string count;
  std::string matrix;
  std::string pa_matrix;
  std::string format;
  std::string id;
  bool sorted {false};
  bool lz4 {false};
  bool lz4_in {false};
  bool no_count {false};
  std::string output;

  std::string display()
  {
    std::stringstream ss;
    ss << this->global_display();
    RECORD(ss, count);
    RECORD(ss, matrix);
    RECORD(ss, pa_matrix);
    RECORD(ss, format);
    RECORD(ss, id);
    RECORD(ss, sorted);
    RECORD(ss, lz4);
    RECORD(ss, lz4_in);
    RECORD(ss, output);
    return ss.str();
  }
};

using agg_options_t = std::shared_ptr<struct agg_options>;

};
