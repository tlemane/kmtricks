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

#include <iostream>

// int
#include <kmtricks/config.hpp>
#include <kmtricks/utils.hpp>
#include <limits>

namespace km
{
inline void main_infos(std::ostream& stream)
{
  stream << "kmtricks " << PROJECT_VER << "\n\n";
  stream << "- HOST -" << "\n";
  stream << "build host: " << HOST_SYSTEM << "\n";
  stream << "run host: " << get_uname_sr() << "\n";
  stream << "- BUILD -" << "\n";
  stream << "c compiler: " << COMPILER_C << "\n";
  stream << "cxx compiler: " << COMPILER_CXX << "\n";
  stream << "conda: " << CONDA_BUILD << "\n";
  stream << "static: " << STATIC_BUILD << "\n";
  stream << "native: " << NATIVE_BUILD << "\n";
  stream << "modules: " << MODULES_BUILD << "\n";
  stream << "socks: " << SOCKS_ON << "\n";
  stream << "howde: " << HOWDE_BUILD << "\n";
  stream << "dev: " << DEV_BUILD << "\n";
  stream << "kmer: " << KMER_LIST_STR << "\n";
  stream << "max_c: " << std::to_string(std::numeric_limits<selectC<DMAX_C>::type>::max()) << "\n";
  stream << "\n";
  stream << "- GIT SHA1 / VERSION -" << "\n";
  stream << "kmtricks: " << GIT_SHA1 << "\n";
  stream << "sdsl: " << SDSL_SHA1 << "\n";
  stream << "bcli: " << BCLI_SHA1 << "\n";
  stream << "fmt: " << FMT_SHA1 << "\n";
  stream << "kff: " << KFF_SHA1 << "\n";
  stream << "lz4: " << LZ4_SHA1 << "\n";
  stream << "spdlog: " << SPDLOG_SHA1 << "\n";
  stream << "xxhash: " << XXHASH_SHA1 << "\n";
  stream << "gtest: " << GTEST_SHA1 << "\n";
  stream << "croaring: " << CROAR_SHA1 << "\n";
  stream << "robin-hood-hasing: " << ROBIN_SHA1 << "\n";
  stream << "turbop: " << TURBOP_SHA1 << "\n";
  stream << "cfrcat: " << CFR_SHA1 << "\n";
  stream << "indicators: " << IND_SHA1 << "\n\n";
  stream << "Contact: " << CONTACT << "\n";
  stream << std::flush;
}

};  // namespace kmdiff