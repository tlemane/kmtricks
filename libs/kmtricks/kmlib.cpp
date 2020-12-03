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
#include "kmlib.hpp"

namespace km
{
template class Code<uint8_t>;
template class Kmer<uint8_t>;
template class Superk<uint8_t>;
template class Hasher<uint8_t>;
template class Validator<uint8_t>;
template class SuperkReader<uint8_t>;

template class Code<uint16_t>;
template class Kmer<uint16_t>;
template class Superk<uint16_t>;
template class Hasher<uint16_t>;
template class Validator<uint16_t>;
template class SuperkReader<uint16_t>;

template class Code<uint32_t>;
template class Kmer<uint32_t>;
template class Superk<uint32_t>;
template class Hasher<uint32_t>;
template class Validator<uint32_t>;
template class Minimizer<uint32_t>;
template class SuperkReader<uint32_t>;

template class Code<uint64_t>;
template class Kmer<uint64_t>;
template class Superk<uint64_t>;
template class Hasher<uint64_t>;
template class Validator<uint64_t>;
template class Minimizer<uint64_t>;
template class SuperkReader<uint64_t>;

template class Merger<uint8_t, uint8_t>;
template class Merger<uint8_t, uint16_t>;
template class Merger<uint8_t, uint32_t>;
template class Merger<uint8_t, uint64_t>;

template class Merger<uint16_t, uint8_t>;
template class Merger<uint16_t, uint16_t>;
template class Merger<uint16_t, uint32_t>;
template class Merger<uint16_t, uint64_t>;

template class Merger<uint32_t, uint8_t>;
template class Merger<uint32_t, uint16_t>;
template class Merger<uint32_t, uint32_t>;
template class Merger<uint32_t, uint64_t>;

template class Merger<uint64_t, uint8_t>;
template class Merger<uint64_t, uint16_t>;
template class Merger<uint64_t, uint32_t>;
template class Merger<uint64_t, uint64_t>;

#ifdef __SIZEOF_INT128__
template class Code<__uint128_t>;
template class Kmer<__uint128_t>;
template class Superk<__uint128_t>;
template class Hasher<__uint128_t>;
template class Validator<__uint128_t>;
template class SuperkReader<__uint128_t>;

template class Merger<uint8_t, __uint128_t>;
template class Merger<uint16_t, __uint128_t>;
template class Merger<uint32_t, __uint128_t>;
template class Merger<uint64_t, __uint128_t>;
template class Merger<__uint128_t, uint8_t>;
template class Merger<__uint128_t, uint16_t>;
template class Merger<__uint128_t, uint32_t>;
template class Merger<__uint128_t, uint64_t>;
template class Merger<__uint128_t, __uint128_t>;
#endif
}