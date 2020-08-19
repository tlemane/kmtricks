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
#include <gatb/gatb_core.hpp>
#include <cmath>
#include <kmtricks/merger.hpp>
#include <kmtricks/bitmatrix.hpp>
#include "config.hpp"

using namespace std;
using namespace km;
typedef unsigned char uchar;

class KmMerge : public gatb::core::tools::misc::impl::Tool
{
public:
  KmMerge();
  void        execute() override;

private:
  void        parse_args();
  void        merge_to_ascii();
  void        merge_to_bin();
  void        merge_to_pa_matrix();
  void        merge_to_bf_pa();
  void        transpose();

private:
  Merger<kmtype_t, cntype_t>  *_m;
  Env*                  e;
  uint                  _min_a;
  uint                  _min_r;
  uint                  _id;
  uint64_t              _lower_hash;
  uint64_t              _upper_hash;
  string                _run_dir;
  string                _fofpath;
  uint8_t               _mode;
  vector<tuple<uint64_t, uint64_t>> _hash_windows;
};

