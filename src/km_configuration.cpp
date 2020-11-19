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

#include <csignal>
#include <fmt/format.h>
#include <kmtricks/logging.hpp>
#include "km_configuration.hpp"
#include "config.hpp"

km::log_config km::LOG_CONFIG;

Kmtricks::Kmtricks()
  : Tool("kmtricks"),
  _min_hash(0) 
{
  setParser(new OptionsParser("kmtricks: build runtime environment"));

  getParser()->push_back(new OptionOneParam(STR_URI_FILE,
                                            "fof that contains path of read files, one per line", true));
  getParser()->push_back(new OptionOneParam(STR_RUN_DIR,
                                            "directory to write tmp and output files", true));
  getParser()->push_back(new OptionOneParam(STR_KMER_SIZE,
                                            "size of a kmer", false, "31"));
  getParser()->push_back(new OptionOneParam(STR_KMER_ABUNDANCE_MIN,
                                            "min abundance threshold for solid kmers", false, "2"));
  getParser()->push_back(new OptionOneParam(STR_KMER_ABUNDANCE_MAX,
                                            "max abundance threshold for solid kmers", false, "max"));
  getParser()->push_back(new OptionOneParam(STR_MAX_MEMORY,
                                            "max memory available in megabytes", false, "8000"));
  getParser()->push_back(new OptionOneParam(STR_NB_CORES, "not used, needed by gatb args parser", false, "1"), 0UL, false);

  IOptionsParser *mParser = new OptionsParser("advanced performance tweaks");
  mParser->push_back(new OptionOneParam(STR_MINIMIZER_TYPE,
                                        "minimizer type (0=lexi, 1=freq)", false, "0"));
  mParser->push_back(new OptionOneParam(STR_MINIMIZER_SIZE,
                                        "size of a minimizer", false, "10"));
  mParser->push_back(new OptionOneParam(STR_REPARTITION_TYPE,
                                        "minimizer repartition (0=unordered, 1=ordered)", false, "0"));
  mParser->push_back(new OptionOneParam(STR_NB_PARTS, "number of partitions (0=auto)", false, "0"));
  getParser()->push_back(mParser);

  IOptionsParser *oParser = new OptionsParser("hash mode configuration, if you want to use kmtricks in hash mode");

  oParser->push_back(new OptionOneParam(STR_HASHER,
                                        "hash function: sabuhash, xor", false, "xor"));
  oParser->push_back(new OptionOneParam(STR_MAX_HASH,
                                        "max hash value ( 0 < hash < max(int64) )", false, "1000000000"));
  getParser()->push_back(oParser);

  km::LOG_CONFIG.show_labels=true;
  km::LOG_CONFIG.level=km::DEBUG;
}

void Kmtricks::parse_args()
{
  _fof_path         = getInput()->getStr(STR_URI_FILE);
  _max_memory       = getInput()->getInt(STR_MAX_MEMORY);
  _k_size           = getInput()->getInt(STR_KMER_SIZE);
  _a_min            = getInput()->getInt(STR_KMER_ABUNDANCE_MIN);
  if (getInput()->getStr(STR_KMER_ABUNDANCE_MAX) == "max")
    _a_max = maxc.at(sizeof(cntype_t));
  else
    _a_max            = getInput()->getInt(STR_KMER_ABUNDANCE_MAX);
  _dir              = getInput()->getStr(STR_RUN_DIR);
  _nb_partitions    = getInput()->getInt(STR_NB_PARTS);

  _max_hash = getInput()->getInt(STR_MAX_HASH);
  _hasher   = getInput()->getStr(STR_HASHER);

  getInput()->add(1, STR_MAX_DISK, "0");
  getInput()->add(1, STR_STORAGE_TYPE, "file");
  getInput()->add(1, STR_SOLIDITY_KIND, "sum");
}

void Kmtricks::init()
{
  e = new Env(_dir, "");
  e->build();
  _f_log = ofstream(e->LOG_CMD, ios::out);

  const size_t span = KMER_DEFAULT_SPAN;

  fof_t fof = parse_km_fof(_fof_path);
  string in = all_files(fof);
  
  IBank* bank = Bank::open(in);
  LOCAL(bank);
  Storage *_config_storage = StorageFactory(STORAGE_FILE).create(e->STORE_CONFIG, true, false);
  LOCAL(_config_storage);
  ConfigurationAlgorithm<span> configA (bank, getInput());
  configA.execute();
  Configuration _config = configA.getConfiguration();
  if (_nb_partitions != 0)
    _config._nb_partitions = _nb_partitions;
  _config.save(_config_storage->getGroup("config"));
  _nb_partitions = _config._nb_partitions;
  e->build_p(_nb_partitions);

  auto window_size = NMOD8((uint64_t)ceil((double)_max_hash / (double)_nb_partitions));

  ofstream hw(e->HASHW_MAP);
  hw.write((char*)&_nb_partitions, sizeof(uint32_t));
  uint64_t lb, ub;
  for ( int i = 0; i < _nb_partitions; i++ )
  {
    lb = i * window_size;
    ub = ((i + 1) * window_size) - 1;
    _hash_windows.emplace_back(make_tuple(lb, ub));
    hw.write((char *) &lb, sizeof(uint64_t));
    hw.write((char *) &ub, sizeof(uint64_t));
  }
  hw.write((char*) &_max_hash, sizeof(uint64_t));
  uint32_t minimsize = _config._minim_size;
  hw.write((char*) &minimsize, sizeof(uint32_t));
  hw.close();

  ofstream log_file(e->DIR+"/config.log", ios::out);
  km::LOG(km::INFO, log_file) << "Fof path:          " << _fof_path;
  km::LOG(km::INFO, log_file) << "Kmer size:         " << _k_size;
  km::LOG(km::INFO, log_file) << "Abundance min/max: " << _a_min << " " << _a_max;
  km::LOG(km::INFO, log_file) << "Nb partitions:     " << _nb_partitions;
  km::LOG(km::INFO, log_file) << "Max hash / hasher: " << _max_hash << " " << _hasher;
  km::LOG(km::INFO, log_file) << "Window size:       " << window_size;
  km::LOG(km::INFO, log_file) << "Minimizer size:    " << minimsize;
  log_file.close();
}

void Kmtricks::execute()
{
  parse_args();
  init();
}

int main(int argc, char* argv[])
{
  try 
  {
    Kmtricks().run(argc, argv);
  }
  catch (OptionFailure &e)
  {
    return e.displayErrors(cout);
  }
  catch (Exception &e)
  {
    cerr << e.getMessage() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
