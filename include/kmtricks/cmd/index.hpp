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

struct index_options : km_options
{
  bool howde;
  bool allsome;
  bool determined;
  bool brief;
  bool uncompressed;
  bool rrr;
  bool roar;
  uint64_t bits;
  double cull;
  bool keep;
  size_t lower;
  size_t upper;
  bool cull2;
  double cullsd;

  std::string display()
  {
    std::stringstream ss;
    ss << this->global_display();
    RECORD(ss, dir);
    RECORD(ss, allsome);
    RECORD(ss, determined);
    RECORD(ss, brief);
    RECORD(ss, uncompressed);
    RECORD(ss, rrr);
    RECORD(ss, roar);
    RECORD(ss, bits);
    RECORD(ss, cull);
    RECORD(ss, keep);
    RECORD(ss, lower);
    RECORD(ss, upper);
    RECORD(ss, cull2);
    RECORD(ss, cullsd);
    std::string ret = ss.str(); ret.pop_back(); ret.pop_back();
    return ret;
  }
};

using index_options_t = std::shared_ptr<struct index_options>;

};
