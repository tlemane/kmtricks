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
#include <thread>
// ext
#include <bcli/bcli.hpp>
// int
#include <kmtricks/cmd/cmd_common.hpp>

namespace km
{

using cli_t = std::shared_ptr<bc::Parser<1>>;

inline void add_common(bc::cmd_t cmd, km_options_t options)
{
  cmd->add_group("common", "");
  cmd->add_param("-t/--threads", "number of threads.")
    ->def(std::to_string(std::thread::hardware_concurrency()))
    ->meta("INT")
    ->setter(options->nb_threads)
    ->checker(bc::check::is_number);
  cmd->add_param("-h/--help", "show this message and exit.")
    ->as_flag()
    ->action(bc::Action::ShowHelp);
  cmd->add_param("--version", "show version and exit.")
    ->as_flag()
    ->action(bc::Action::ShowVersion);
  cmd->add_param("-v/--verbose", "verbosity level [debug|info|warning|error].")
    ->meta("STR")
    ->def("info")
    ->checker(bc::check::f::in("debug|info|warning|error"))
    ->setter(options->verbosity);
}

};  // namespace kmdiff