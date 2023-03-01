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

#include <string>
#include <vector>

#include <kmtricks/cli/cli_common.hpp>
#include <kmtricks/cmd/cmd_common.hpp>

namespace km {

struct combine_options : km_options
{
  std::vector<std::string> runs;
  std::string output;
  bool cpr;

  std::string display()
  {
    std::stringstream ss;
    ss << this->global_display();
    RECORD(ss, output);
    RECORD(ss, cpr);
    std::string ret = ss.str(); ret.pop_back(); ret.pop_back();
    return ret;
  }
};

using combine_options_t = std::shared_ptr<struct combine_options>;

};
