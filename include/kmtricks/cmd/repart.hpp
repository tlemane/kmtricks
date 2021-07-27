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

struct repart_options : km_options
{
  std::string fof;
  uint32_t kmer_size;

  uint32_t minim_type;
  uint32_t minim_size;
  uint32_t repart_type;
  uint32_t nb_parts;
  uint64_t bloom_size;

  std::string display()
  {
    std::stringstream ss;
    ss << this->global_display();
    RECORD(ss, fof);
    RECORD(ss, kmer_size);
    RECORD(ss, minim_size);
    RECORD(ss, minim_type);
    RECORD(ss, repart_type);
    RECORD(ss, nb_parts);
    std::string ret = ss.str(); ret.pop_back(); ret.pop_back();
    return ret;
  }
};

using repart_options_t = std::shared_ptr<struct repart_options>;

};