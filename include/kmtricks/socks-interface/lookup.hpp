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
#include <kmtricks/cmd/cmd_common.hpp>
#include <memory>
#include <thread>

#include <spdlog/spdlog.h>
#include <kmtricks/cmd/cmd_common.hpp>
#include <kmtricks/cli/cli_common.hpp>
#include <kmtricks/cmd.hpp>
#include <kmtricks/socks-interface/socks_utils.hpp>

namespace km {

struct lookup_options : km_options
{
  std::string output;
  double threshold;
  std::string out_type;
  std::string query;
  int z;

  std::string display()
  {
    std::stringstream ss;
    ss << this->global_display();
    RECORD(ss, output);
    RECORD(ss, out_type);
    RECORD(ss, query);
    RECORD(ss, threshold);
    RECORD(ss, z);
    return ss.str();
  }
};

using lookup_options_t = std::shared_ptr<struct lookup_options>;

inline km_options_t lookup_cli(std::shared_ptr<bc::Parser<1>> cli, lookup_options_t options)
{
  bc::cmd_t look_cmd = cli->add_command("lookup-kmer", "Lookup k-mers.");

  auto pset = [options](const std::string& v) {
    options->query = v;
  };

  look_cmd->set_positionals(1, "", "");
  look_cmd->set_positionals_help("<query>", "A query file in fasta format");
  look_cmd->positionals_setter(pset);

  look_cmd->add_param("-i/--index-dir", "Index directory.")
    ->meta("DIR")
    ->checker(bc::check::is_dir)
    ->setter(options->dir);

  look_cmd->add_param("-t/--threshold", "Threshold.")
    ->meta("FLOAT")
    ->def("0.7")
    ->checker(bc::check::f::range(0.0, 1.0))
    ->setter(options->threshold);

  look_cmd->add_param("-o/--output-type", "output type. [vector|list]")
    ->meta("STR")
    ->def("vector")
    ->checker(bc::check::f::in("vector|list"))
    ->setter(options->out_type);

  add_common(look_cmd, options);
  return options;
}

template<size_t MAX_K>
struct main_lookup
{
  void operator()(km_options_t options)
  {
    lookup_options_t opt = std::static_pointer_cast<struct lookup_options>(options);
    query_options_t query_opt = std::make_shared<struct query_options>(query_options{});

    KmDir::get().init(opt->dir, "", false);
    std::string tmp = fmt::format("{}/tmp_query_res", KmDir::get().m_index_storage);
    query_opt->dir = opt->dir;
    query_opt->query = opt->query;
    query_opt->threshold = opt->threshold;
    query_opt->output = tmp;

    KmDir::get().init(opt->dir, "", false);
    main_query<MAX_K>()(query_opt);

    std::vector<std::string> query_idx;
    {
      std::ifstream in(opt->query, std::ios::in); check_fstream_good(opt->query, in);
      for (std::string line; std::getline(in, line);)
      {
        if (line.find(">") != std::string::npos)
          query_idx.push_back(line.substr(1));
      }
    }

    if (opt->out_type == "vector")
    {
      format_result_vector(tmp, std::cout, query_idx, KmDir::get().m_fof);
    }
    else
    {
      format_result_list(tmp, std::cout, query_idx, KmDir::get().m_fof);
    }
  }
};

};