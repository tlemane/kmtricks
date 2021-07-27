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

#define KMTRICKS_PUBLIC
#include <kmtricks/kmer.hpp>
#include <kmtricks/minimizer.hpp>
#include <kmtricks/repartition.hpp>
#include <kmtricks/histogram.hpp>
#include <kmtricks/bitmatrix.hpp>
#include <kmtricks/hash.hpp>
#include <kmtricks/timer.hpp>
#include <kmtricks/itask.hpp>
#include <kmtricks/hash.hpp>
#include <kmtricks/loop_executor.hpp>
#include <kmtricks/task_pool.hpp>

#ifdef WITH_KM_IO // requires lz4 and turpob
#include <kmtricks/io.hpp>
#include <kmtricks/merge.hpp>
#endif

#ifdef WITH_KM_HOWDE // requires sdsl and km_howdesbt/bloom_filter_file.h
#include <kmtricks/howde_utils.hpp>
#endif