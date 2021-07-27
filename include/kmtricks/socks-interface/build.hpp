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
#include <kmtricks/config.hpp>
#include <kmtricks/cmd/cmd_common.hpp>
#include <kmtricks/cli/cli_common.hpp>
#include <kmtricks/cmd/all.hpp>
#include <kmtricks/cmd/index.hpp>
#include <kmtricks/cmd.hpp>
#include <kmtricks/loop_executor.hpp>

namespace km {

struct build_options : km_options
{
  std::string input;
  std::string output;
  size_t kmer_size {0};
  size_t ab_min {0};
  size_t bloom_size {0};
  size_t nb_parts {0};
  std::vector<std::string> pos;
  std::string display()
  {
    std::stringstream ss;
    ss << this->global_display();
    return ss.str();
  }
};

using build_options_t = std::shared_ptr<struct build_options>;

km_options_t build_cli(std::shared_ptr<bc::Parser<1>> cli, build_options_t options)
{
  bc::cmd_t build_cmd = cli->add_command("build", "Build index.");

  auto pset = [options](const std::string& v) {
    options->pos.push_back(v);
  };

  build_cmd->set_positionals(2, "", "");
  build_cmd->set_positionals_help("<input> <output>",
                                  "<input> : kmtricks fof\n  <output> : output directory");
  build_cmd->positionals_setter(pset);


  build_cmd->add_param("-k/--kmer-size", "size of k-mers.")
    ->meta("INT")
    ->def("31")
    ->checker(bc::check::is_number)
    ->setter(options->kmer_size);

  build_cmd->add_param("-m/--abundance-min", "Min abundance for solid k-mers.")
    ->meta("INT")
    ->def("1")
    ->checker(bc::check::is_number)
    ->setter(options->ab_min);

  build_cmd->add_param("-b/--bloom-size", "size of k-mers.")
    ->meta("INT")
    ->def("100000000")
    ->checker(bc::check::is_number)
    ->setter(options->bloom_size);

  build_cmd->add_param("--nb-partitions", "number of partitions (0=auto).")
    ->meta("INT")
    ->def("0")
    ->checker(bc::check::is_number)
    ->setter(options->nb_parts);

  add_common(build_cmd, options);
  return options;
}

template<size_t MAX_K>
struct main_build
{
  void operator()(km_options_t options)
  {
    build_options_t opt = std::static_pointer_cast<struct build_options>(options);
    all_options_t all_opt = std::make_shared<struct all_options>(all_options{});
    index_options_t index_opt = std::make_shared<struct index_options>(index_options{});
    opt->input = opt->pos[0]; opt->output = opt->pos[1];
    all_opt->dir = opt->output;
    all_opt->fof = opt->input;
    all_opt->kmer_size = opt->kmer_size;
    all_opt->c_ab_min = opt->ab_min;
    all_opt->bloom_size = opt->bloom_size;
    all_opt->format = FORMAT::BIN;
    all_opt->mode = MODE::BFT;
    all_opt->count_format = COUNT_FORMAT::HASH;
    all_opt->out_format = OUT_FORMAT::HOWDE;
    all_opt->minim_size = 10;
    all_opt->restrict_to = 1.0;
    all_opt->nb_threads = opt->nb_threads;
    all_opt->nb_parts = opt->nb_parts;
    all_opt->skip_merge = true;
    main_all<MAX_K>()(all_opt);

    index_opt->dir = opt->output;
    index_opt->howde = true;
    index_opt->bits = opt->bloom_size * 0.1;
    main_index<MAX_K>()(index_opt);
  }
};


};