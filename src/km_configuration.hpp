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
#include <libgen.h>
#include <gatb/gatb_core.hpp>
#include <gatb/kmer/impl/ConfigurationAlgorithm.hpp>
#include "config.hpp"
#include "kmtricks/utilities.hpp"

class Kmtricks : public Tool
{
public:
    Kmtricks();
private:
  void parse_args();
  void init();
  void execute() override;

  IteratorListener* _progress;
  vector<string>    _bank_paths;
  ofstream          _f_log;
  
  string            _fof_path;
  Env*              e;
  size_t            _k_size;
  uint              _a_min, _a_max;
  uint              _nb_cores, _max_memory;
  uint              _nb_partitions;
  uint              _nb_procs;

  string            _dir;

  string            _hasher;
  uint64_t          _min_hash, _max_hash;


  vector<tuple<uint64_t, uint64_t>> _hash_windows;
};

