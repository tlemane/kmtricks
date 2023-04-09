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

// ext
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

// int
#include <kmtricks/cmd.hpp>
#include <kmtricks/cli.hpp>
#include <kmtricks/signals.hpp>

#include <kmtricks/loop_executor.hpp>

#include <gatb/gatb_core.hpp>

int main(int argc, char* argv[])
{
  using namespace km;

  SignalHandler::get().init();
  kmtricksCli cli(PROJECT_NAME, PROJECT_DESC, PROJECT_VER, "");

  auto [cmd, options] = cli.parse(argc, argv);

  set_verbosity_level(options->verbosity);
  auto cerr_logger = spdlog::stderr_color_mt("kmtricks");
  cerr_logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
  spdlog::set_default_logger(cerr_logger);

  size_t kmer_size;
  if (cmd != COMMAND::ALL && cmd != COMMAND::REPART && cmd != COMMAND::INFOS)
  {
    KmDir::get().init(options->dir, "", false);
    Storage* config_storage = StorageFactory(STORAGE_FILE).load(KmDir::get().m_config_storage);
    LOCAL(config_storage);
    Configuration config = Configuration();
    config.load(config_storage->getGroup("gatb"));
    kmer_size = config._kmerSize;
  }
  else if (cmd == COMMAND::REPART)
  {
    kmer_size = std::static_pointer_cast<struct repart_options>(options)->kmer_size;
  }
  else if (cmd == COMMAND::ALL)
  {
    kmer_size = std::static_pointer_cast<struct all_options>(options)->kmer_size;
  }

  try
  {
    if (cmd == COMMAND::ALL)
    {
      const_loop_executor<0, KMER_N>::exec<main_all>(kmer_size, options);
    }
    else if (cmd == COMMAND::REPART)
    {
      const_loop_executor<0, KMER_N>::exec<main_repart>(kmer_size, options);
    }
    else if (cmd == COMMAND::SUPERK)
    {
      const_loop_executor<0, KMER_N>::exec<main_superk>(kmer_size, options);
    }
    else if (cmd == COMMAND::COUNT)
    {
      const_loop_executor<0, KMER_N>::exec<main_count>(kmer_size, options);
    }
    else if (cmd == COMMAND::MERGE)
    {
      const_loop_executor<0, KMER_N>::exec<main_merge>(kmer_size, options);
    }
    else if (cmd == COMMAND::FORMAT)
    {
      const_loop_executor<0, KMER_N>::exec<main_format>(kmer_size, options);
    }
    else if (cmd == COMMAND::DUMP)
    {
      const_loop_executor<0, KMER_N>::exec<main_dump>(kmer_size, options);
    }
    else if (cmd == COMMAND::AGGREGATE)
    {
      const_loop_executor<0, KMER_N>::exec<main_agg>(kmer_size, options);
    }
    else if (cmd == COMMAND::FILTER)
    {
      const_loop_executor<0, KMER_N>::exec<main_filter>(kmer_size, options);
    }
    else if (cmd == COMMAND::COMBINE)
    {
      const_loop_executor<0, KMER_N>::exec<main_combine>(kmer_size, options);
    }
#ifdef WITH_HOWDE
    else if (cmd == COMMAND::INDEX)
    {
      const_loop_executor<0, KMER_N>::exec<main_index>(kmer_size, options);
    }
    else if (cmd == COMMAND::QUERY)
    {
      const_loop_executor<0, KMER_N>::exec<main_query>(kmer_size, options);
    }
#endif
    else if (cmd == COMMAND::INFOS)
    {
      main_infos(std::cerr);
    }

  }
  catch (const km_exception& e)
  {
    spdlog::error("{} - {}", e.get_name(), e.get_msg());
    exit(EXIT_FAILURE);
  }
  catch (const Exception& e)
  {
    spdlog::error("GATB ERROR: {}", e.getMessage());
    exit(EXIT_FAILURE);
  }
  catch (const std::exception& e)
  {
    spdlog::error(e.what());
    exit(EXIT_FAILURE);
  }

  return 0;
}
