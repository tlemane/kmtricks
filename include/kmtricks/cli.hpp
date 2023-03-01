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

// std
#include <memory>
#include <string>

// ext
#include <bcli/bcli.hpp>

// int
#include <kmtricks/config.hpp>
#include <kmtricks/cli/cli_common.hpp>
#include <kmtricks/cli/all.hpp>
#include <kmtricks/cli/repart.hpp>
#include <kmtricks/cli/superk.hpp>
#include <kmtricks/cli/count.hpp>
#include <kmtricks/cli/merge.hpp>
#include <kmtricks/cli/format.hpp>
#include <kmtricks/cli/dump.hpp>
#include <kmtricks/cli/infos.hpp>
#include <kmtricks/cli/aggregate.hpp>
#include <kmtricks/cli/filter.hpp>
#include <kmtricks/cli/index.hpp>
#include <kmtricks/cli/query.hpp>
#include <kmtricks/cli/combine.hpp>

namespace km
{

constexpr int KL[KMER_N] = {KMER_LIST};

class kmtricksCli
{
public:
  kmtricksCli(
      const std::string& name,
      const std::string& desc,
      const std::string& version,
      const std::string& authors);

  std::tuple<COMMAND, km_options_t> parse(int argc, char* argv[]);

private:
  cli_t cli {nullptr};
  all_options_t all_opt {nullptr};
  repart_options_t repart_opt {nullptr};
  superk_options_t superk_opt {nullptr};
  count_options_t count_opt {nullptr};
  merge_options_t merge_opt {nullptr};
  format_options_t format_opt {nullptr};
  dump_options_t dump_opt {nullptr};
  agg_options_t agg_opt {nullptr};
  filter_options_t filter_opt {nullptr};
  index_options_t index_opt {nullptr};
  query_options_t query_opt {nullptr};
  combine_options_t combine_opt {nullptr};
};

};  // namespace km
