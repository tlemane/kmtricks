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
#include <sdsl/bit_vectors.hpp>
#include <HowDeSBT/bloom_filter_file.h>
#include "config.hpp"

#define NBYTE(bits) (((bits) >> 3) + ((bits) % 8 != 0))
#define round_up_16(b)  ((((std::uint64_t) (b))+15)&(~15))
typedef sdsl::bit_vector bitvector;

class KmConvert : public Tool
{
public:
  KmConvert ();

private:
  void execute();
  void parse_args();
  void init();

private:
  Env*    _e;
  string  _run_dir;
  string  _split_str;
  string  _hm_path;
  bool    _howde;
  bool     _sdsl;
  uint64_t _vlen;
  uint64_t    _filter_size;
  uint    _nb_files;
  uint    _nb_parts;
  uint64_t _win_size;
  uint32_t _kmer_size;
  string _fof;
  string _sync;
  vector<string>   _f_names;
  vector<ifstream> _matrices;
  vector<tuple<uint64_t, uint64_t>> _hash_windows;
  map<string, uint> _fof_pos;
  vector<uint>      _pos;


};