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

#include <kmtricks/config.hpp>
#include <kmtricks/signals.hpp>
#include <kmtricks/socks-interface/build.hpp>
#include <kmtricks/socks-interface/lookup.hpp>
#include <kmtricks/kmdir.hpp>
#include <kmtricks/exceptions.hpp>
#include <kmtricks/loop_executor.hpp>
#include <gatb/gatb_core.hpp>

namespace km {

class socksCli
{
public:
  socksCli(const std::string& name, const std::string& desc, const std::string& version,
           const std::string& authors)
  {
    cli = std::make_shared<bc::Parser<1>>(bc::Parser<1>(name, desc, version, authors));
    build_opt = std::make_shared<struct build_options>(build_options{});
    look_opt = std::make_shared<struct lookup_options>(lookup_options{});
    build_cli(cli, build_opt);
    lookup_cli(cli, look_opt);
  }

  std::tuple<COMMAND, km_options_t> parse(int argc, char* argv[])
  {
    try
    {
      (*cli).parse(argc, argv);
    }
    catch (const bc::ex::BCliError& e)
    {
      bc::utils::exit_bcli(e);
      exit(EXIT_FAILURE);
    }

    if (cli->is("build"))
      return std::make_tuple(COMMAND::SOCKS_BUILD, build_opt);
    else if (cli->is("lookup-kmer"))
      return std::make_tuple(COMMAND::SOCKS_LOOKUP, look_opt);
    else
      return std::make_tuple(COMMAND::INFOS, std::make_shared<struct km_options>(km_options{}));
  }

private:
  cli_t cli {nullptr};
  build_options_t build_opt {nullptr};
  lookup_options_t look_opt {nullptr};
};

};

int main(int argc, char* argv[])
{
  using namespace km;

  SignalHandler::get().init();
  socksCli cli("kmtricks-socks", "kmtricks socks interface", PROJECT_VER, "");

  auto [cmd, options] = cli.parse(argc, argv);

  set_verbosity_level(options->verbosity);
  auto cerr_logger = spdlog::stderr_color_mt("kmtricks");
  cerr_logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
  spdlog::set_default_logger(cerr_logger);

  size_t kmer_size;

  try
  {
    if (cmd == COMMAND::SOCKS_BUILD)
    {
      kmer_size = std::static_pointer_cast<struct build_options>(options)->kmer_size;
      const_loop_executor<0, KMER_N>::exec<main_build>(kmer_size, options);
    }
    else if (cmd == COMMAND::SOCKS_LOOKUP)
    {
      KmDir::get().init(options->dir, "", false);
      Storage* config_storage = StorageFactory(STORAGE_FILE).load(KmDir::get().m_config_storage);
      LOCAL(config_storage);
      Configuration config = Configuration();
      config.load(config_storage->getGroup("gatb"));
      kmer_size = config._kmerSize;
      const_loop_executor<0, KMER_N>::exec<main_lookup>(kmer_size, options);
    }

  }
  catch (const km_exception& e)
  {
    spdlog::error("{} - {}", e.get_name(), e.get_msg());
  }
  catch (const Exception& e)
  {
    spdlog::error("GATB ERROR: {}", e.getMessage());
  }
  catch (const std::exception& e)
  {
    spdlog::error(e.what());
  }

  return 0;
}