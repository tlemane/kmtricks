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
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#define RECORD(ss, var) ss << #var << "=" << var << ", "

namespace km
{

enum class COMMAND
{
  ALL,
  REPART,
  SUPERK,
  COUNT,
  MERGE,
  FORMAT,
  DUMP,
  AGGREGATE,
  FILTER,
  INDEX,
  QUERY,
  INFOS,
  SOCKS_BUILD,
  SOCKS_LOOKUP,
  COMBINE,
  UNKNOWN
};

inline COMMAND str_to_cmd(const std::string& s)
{
  if (s == "repart")
    return COMMAND::REPART;
  else if (s == "superk")
    return COMMAND::SUPERK;
  else if (s == "count")
    return COMMAND::COUNT;
  else if (s == "merge")
    return COMMAND::MERGE;
  else if (s == "format")
    return COMMAND::FORMAT;
  else if (s == "dump")
    return COMMAND::DUMP;
  else if (s == "aggregate")
    return COMMAND::AGGREGATE;
  else if (s == "filter")
    return COMMAND::FILTER;
  else if (s == "index")
    return COMMAND::INDEX;
  else if (s == "query")
    return COMMAND::QUERY;
  else if (s == "build")
    return COMMAND::SOCKS_BUILD;
  else if (s == "lookup-kmer")
    return COMMAND::SOCKS_LOOKUP;
  else if (s == "combine")
    return COMMAND::COMBINE;
  else
    return COMMAND::ALL;
}

inline std::string cmd_to_str(COMMAND cmd)
{
  if (cmd == COMMAND::REPART)
    return "repart";
  else if (cmd == COMMAND::SUPERK)
    return "superk";
  else if (cmd == COMMAND::COUNT)
    return "count";
  else if (cmd == COMMAND::MERGE)
    return "merge";
  else if (cmd == COMMAND::FORMAT)
    return "format";
  else if (cmd == COMMAND::DUMP)
    return "dump";
  else if (cmd == COMMAND::AGGREGATE)
    return "aggregate";
  else if (cmd == COMMAND::FILTER)
    return "filter";
  else if (cmd == COMMAND::INDEX)
    return "index";
  else if (cmd == COMMAND::QUERY)
    return "query";
  else if (cmd == COMMAND::SOCKS_BUILD)
    return "socks-build";
  else if (cmd == COMMAND::SOCKS_LOOKUP)
    return "socks-lookup";
  else if (cmd == COMMAND::COMBINE)
    return "combine";
  else
    return "all";
}

enum class MODE
{
  COUNT,
  TEXT,
  BIN,
  PA,
  BF,
  BFT,
  BFC,
  UNKNOWN,
};

inline MODE str_to_mode(const std::string& s)
{
  if (s == "bin")
    return MODE::BIN;
  else if (s == "text")
    return MODE::TEXT;
  else if (s == "count")
    return MODE::COUNT;
  else if (s == "pa")
    return MODE::PA;
  else if (s == "bf")
    return MODE::BF;
  else if (s == "bft")
    return MODE::BFT;
  else if (s == "bfc")
    return MODE::BFC;
  else
    return MODE::UNKNOWN;
}

inline std::string mode_to_str(MODE mode)
{
  if (mode == MODE::BIN)
    return "bin";
  else if (mode == MODE::TEXT)
    return "text";
  else if (mode == MODE::COUNT)
    return "count";
  else if (mode == MODE::PA)
    return "pa";
  else if (mode == MODE::BF)
    return "bf";
  else if (mode == MODE::BFT)
    return "bft";
  else if (mode == MODE::BFC)
    return "bfc";
  else
    return "unknown";
}

enum class HASHER
{
  XOR,
  XXHASH,
  SABUHASH,
  UNKNOWN
};

inline HASHER str_to_hasher(const std::string& s)
{
  if (s == "XOR")
    return HASHER::XOR;
  else if (s == "XXHASH")
    return HASHER::XXHASH;
  else if (s == "SABUHASH")
    return HASHER::SABUHASH;
  else
    return HASHER::UNKNOWN;
}

enum class OUT_FORMAT
{
  RAW,
  HOWDE,
  SDSL,
  UNKNOWN
};

inline OUT_FORMAT str_to_format(const std::string& s)
{
  if (s == "raw")
    return OUT_FORMAT::RAW;
  else if (s == "howdesbt")
    return OUT_FORMAT::HOWDE;
  else if (s == "sdsl")
    return OUT_FORMAT::SDSL;
  else
    return OUT_FORMAT::UNKNOWN;
}

inline std::string format_to_str(OUT_FORMAT format)
{
  if (format == OUT_FORMAT::RAW)
    return "raw";
  else if (format == OUT_FORMAT::HOWDE)
    return "howdesbt";
  else if (format == OUT_FORMAT::SDSL)
    return "sdsl";
  else
    return "unknown";
}

enum class COUNT_FORMAT
{
  KMER,
  HASH,
  UNKNOWN
};

inline COUNT_FORMAT str_to_cformat(const std::string& s)
{
  if (s == "kmer")
    return COUNT_FORMAT::KMER;
  else if (s == "hash")
    return COUNT_FORMAT::HASH;
  else
    return COUNT_FORMAT::UNKNOWN;
}

inline std::string cformat_to_str(COUNT_FORMAT format)
{
  if (format == COUNT_FORMAT::KMER)
    return "kmer";
  else if (format == COUNT_FORMAT::HASH)
    return "hash";
  else
    return "unknown";
}

enum class FORMAT
{
  BIN,
  TEXT,
  UNKNOWN
};

inline FORMAT str_to_format2(const std::string& s)
{
  if (s == "text")
    return FORMAT::TEXT;
  else if (s == "bin")
    return FORMAT::BIN;
  else
    return FORMAT::UNKNOWN;
}

inline std::string format_to_str2(FORMAT format)
{
  if (format == FORMAT::TEXT)
    return "text";
  else if (format == FORMAT::BIN)
    return "bin";
  else
    return "unknown";
}

struct km_options
{
  std::string verbosity{};
  int nb_threads{1};
  std::string dir;

  std::string global_display()
  {
    std::stringstream ss;
    ss << "Options: ";
    RECORD(ss, dir);
    RECORD(ss, verbosity);
    RECORD(ss, nb_threads);
    return ss.str();
  }
};

using km_options_t = std::shared_ptr<struct km_options>;

};  // namespace kmdiff
