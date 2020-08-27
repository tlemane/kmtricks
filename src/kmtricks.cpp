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
#include "kmtricks.hpp"
#include "synchronizer.hpp"

#if __APPLE__
#include <mach-o/dyld.h>
#endif

void signal_callback(int signum)
{
  cout << "\nInterrupt signal (" << signum << ") received" << endl;
  std::system(KILLALL.c_str());
  cout << "All children are killed." << endl;
  exit(EXIT_FAILURE);
}

Kmtricks::Kmtricks(bool env)
  : Tool("kmtricks"),
  _min_hash(0), _build_runtime(env)
{
  if (!env)
  {
    setParser(new OptionsParser("kmtricks"));
    getParser()->push_back(new OptionOneParam(STR_URI_FILE,
                                              "fof that contains path of read files, one per line", true));
    getParser()->push_back(new OptionOneParam(STR_RUN_DIR,
                                              "directory to write tmp and output files", true));
    getParser()->push_back(new OptionOneParam(STR_KMER_SIZE,
                                              "size of a kmer", false, "31"));
    getParser()->push_back(new OptionOneParam(STR_KMER_ABUNDANCE_MIN,
                                              "min abundance threshold for solid kmers", false, "2"));
    getParser()->push_back(new OptionOneParam(STR_KMER_ABUNDANCE_MAX,
                                              "max abundance threshold for solid kmers", false, "3000000"));
    getParser()->push_back(new OptionOneParam(STR_REC_MIN,
                                              "min recurrence threshold through datasets for solid kmers", false, "2"));
    getParser()->push_back(new OptionOneParam(STR_MAX_MEMORY,
                                              "max memory available in megabytes", false, "8000"));
    getParser()->push_back(new OptionOneParam(STR_MAT_FMT,
                                          "output matrix format: ascii, bin, pa, bf, bf_trp", false, "bin"));
    getParser()->push_back(new OptionOneParam(STR_NB_CORES, "number of cores", false, "8"));
    getParser()->push_back(new OptionOneParam(STR_KEEP_TMP, "keep all files", false, "0"));
    getParser()->push_back(new OptionOneParam(STR_VERBOSE, "verbosity", false, "1"));

    IOptionsParser *cParser = new OptionsParser("kmtricks pipeline control");

    cParser->push_back(new OptionOneParam(STR_UP_TO,
                                              "run until step : part, superk, count, merge", false,
                                              "all"));
    cParser->push_back(new OptionOneParam(STR_ONLY,
                                              "run only step : part, superk, count, merge", false,
                                              "all"));
    getParser()->push_back(cParser);

    IOptionsParser *mParser = new OptionsParser("advanced performance tweaks");

    mParser->push_back(new OptionOneParam(STR_MINIMIZER_TYPE,
                                          "minimizer type (0=lexi, 1=freq)", false, "0"));
    mParser->push_back(new OptionOneParam(STR_MINIMIZER_SIZE,
                                          "size of a minimizer", false, "10"));
    mParser->push_back(new OptionOneParam(STR_REPARTITION_TYPE,
                                          "minimizer repartition (0=unordered, 1=ordered)", false, "0"));
    mParser->push_back(new OptionOneParam(STR_NB_PARTS, "number of partitions", false, "0"));
    getParser()->push_back(mParser);

    IOptionsParser *oParser = new OptionsParser("hash mode configuration, only with -matrix-fmt <bf | bf_trp>");

    oParser->push_back(new OptionOneParam(STR_HASHER,
                                          "hash function: sabuhash, xor", false, "xor"));
    oParser->push_back(new OptionOneParam(STR_MAX_HASH,
                                          "max hash value ( 0 < hash < max(int64) )", false, "1000000000"));
    oParser->push_back(new OptionOneParam(STR_SPLIT,
                                          "split matrix in individual files: sdsl, howde, (only with -matrix-fmt bf_trp)", false, "none"));
    getParser()->push_back(oParser);
  }
  else
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
    getParser()->push_back(new OptionOneParam(STR_NB_CORES, "not used, needed by gatb args parser", false, "1"));

    IOptionsParser *mParser = new OptionsParser("advanced performance tweaks");
    mParser->push_back(new OptionOneParam(STR_MINIMIZER_TYPE,
                                          "minimizer type (0=lexi, 1=freq)", false, "0"));
    mParser->push_back(new OptionOneParam(STR_MINIMIZER_SIZE,
                                          "size of a minimizer", false, "10"));
    mParser->push_back(new OptionOneParam(STR_REPARTITION_TYPE,
                                          "minimizer repartition (0=unordered, 1=ordered)", false, "0"));
    mParser->push_back(new OptionOneParam(STR_NB_PARTS, "number of partitions", false, "0"));
    getParser()->push_back(mParser);

    IOptionsParser *oParser = new OptionsParser("hash mode configuration, if you want to use kmtricks in hash mode");

    oParser->push_back(new OptionOneParam(STR_HASHER,
                                          "hash function: sabuhash, xor", false, "xor"));
    oParser->push_back(new OptionOneParam(STR_MAX_HASH,
                                          "max hash value ( 0 < hash < max(int64) )", false, "1000000000"));
    getParser()->push_back(oParser);
  }
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
  if (!_build_runtime)
  {
    _r_min    = getInput()->getInt(STR_REC_MIN);

    _nb_cores = getInput()->getInt(STR_NB_CORES);
    _mat_str  = getInput()->getStr(STR_MAT_FMT);
    _mat_fmt  = output_format.at(_mat_str);

    _upto = exec_control.at(getInput()->getStr(STR_UP_TO));
    _only = exec_control.at(getInput()->getStr(STR_ONLY));

    _keep_tmp   = getInput()->getInt(STR_KEEP_TMP);
    _str_split  = getInput()->getStr(STR_SPLIT);
    _split      = filter_format.at(_str_split);

    if (_split && (_mat_fmt != 4))
    {
      cout << fmt::format("-split option incompatible with -matrix-fmt {}", _mat_str) << endl;
      cout << getParser()->getHelp() << endl;
      exit(EXIT_FAILURE);
    }
  }

  getInput()->add(1, STR_MAX_DISK, "0");
  getInput()->add(1, STR_STORAGE_TYPE, "file");
  getInput()->add(1, STR_SOLIDITY_KIND, "sum");
}

void Kmtricks::init()
{
  _nb_procs = _nb_cores;

  char *buffer = new char[1024];

#if __APPLE__
  uint32_t size=1024;
  _NSGetExecutablePath(buffer, &size);
#else
  readlink("/proc/self/exe", buffer, 1024);
#endif
  _path_binary = dirname(buffer);
  e = new Env(_dir, _path_binary);
  e->build();
  _f_log = ofstream(e->LOG_CMD, ios::out);

  string line, in;
  ifstream fof(_fof_path);
  if (!fof)
  {
    cout << "Unable to open " << _fof_path << endl; exit(EXIT_FAILURE);
  }
  ofstream copyfof(e->FOF_FILE);
  while( getline(fof, line) )
  {
    _bank_paths.push_back(line);
    in += line + ",";
    copyfof << line << "\n";
  }
  in.pop_back();
  fof.close(); 
  copyfof.close();
  _mode = 0;
  if (_mat_fmt == 3 || _mat_fmt == 4)
    _mode = 1;

  const size_t span = KMER_DEFAULT_SPAN;
  
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
}

void Kmtricks::execute()
{
  if (System::file().doesExistDirectory("./km_backtrace"))
    std::system("rm -rf ./km_backtrace");
  parse_args();
  init();

  signal(SIGINT, signal_callback);

  bool all = (_only == 0);

  if (_build_runtime)
    goto end;

  if ( (_only == 1) | all )
    km_part();
  if (_upto == 1) goto end;

  if ( (_only == 2) | all )
    km_superk();
  if (_upto == 2) goto end;

  if ( (_only == 3) | all )
    km_count();
  if (_upto == 3) goto end;

  if ( (_only == 4) | all )
    km_merger();
  if (_upto == 4) goto end;

  if (_split)
  {
    if ((_only == 5) | all )
      km_output();
    if ( _upto == 5 )
      goto end;
  }
  goto end;
  end:
    _f_log.close();
    delete e;
    cout << endl;
}

void Kmtricks::km_part()
{
  if (!System::file().doesExist(e->STORE + "/partition_storage_gatb/minimRepart.minimRepart"))
  {
    _progress = new ProgressSynchro(
      createIteratorListener(_nb_partitions, "km_minim_repart"), System::thread().newSynchronizer()
    );
    _progress->init();
    string command = fmt::format(PARTITIONER_CMD,
                                 "",
                                 e->PARTITIONER_BIN, // partitioner
                                 _fof_path,          // -file
                                 _k_size,            // -kmer-size
                                 _nb_procs,          // -nb-procs
                                 e->DIR,             // -run-dir
                                 e->LOG_PARTITIONER  // &> partitioner.log
    );
    _f_log << command << endl;
    std::system(command.c_str());
    while (!System::file().doesExist(e->SYNCHRO_P + END_TEMP_P))
    {
      sleep(1);
    }
    _progress->finish();
    delete _progress;
  }
  std::system(fmt::format(RM, e->SYNCHRO_P).c_str());
}

void Kmtricks::km_superk()
{
  uint rmemSk = 300;
  size_t max_cj = min(_nb_procs, (_max_memory/rmemSk));

  _progress = new ProgressSynchro(
    createIteratorListener (_bank_paths.size(), "km_reads_to_superk"), System::thread().newSynchronizer());

  size_t* jobs = new size_t(0);
  FBasedSync _superk_sync = FBasedSync(_bank_paths, e->SYNCHRO_S, END_TEMP_S, jobs, max_cj, _progress, e, 0);
  
  for (size_t i=0; i<_bank_paths.size(); i++)
  {
    string command = fmt::format(SUPERK_CMD,
      "",
      e->SUPERK_BIN,                // superk
      _bank_paths[i],               // -file
      e->DIR,                       // -run-dir
      _k_size,                      // -kmer-size
      1,                            // -nb-cores
      fmt::format(e->LOG_SUPERK, i) // &> superk.log
    );
    _f_log << command << endl;
    std::system(command.c_str());
    ++*jobs;
    _superk_sync.push(_bank_paths[i]);
    _superk_sync.wait_ressources(i);
  }
  _superk_sync.wait_end();
  delete jobs;
  std::system(fmt::format(RM, e->SYNCHRO_S).c_str());
}

void Kmtricks::km_count()
{
  uint64_t nbs = _bank_paths.size()*_nb_partitions;
  _progress = new ProgressSynchro(
    createIteratorListener (nbs, "km_superk_to_kmer_counts"), System::thread().newSynchronizer());

  size_t* jobs = new size_t(0);
  vector<string> count_jobs;
  string jname;
  char cop[256];
  for (size_t i=0; i<_bank_paths.size(); i++)
  {
    strcpy(cop, _bank_paths[i].c_str()); // needed with osx as basename needs a non const char*
    string prefix = basename(cop);
    string name = e->STORE_SUPERK + fmt::format(TEMP_S, prefix);
    PartiInfo<5> tmpinfo(name);
    for (size_t n=0; n<_nb_partitions; n++)
    {
      jname = _bank_paths[i]+"_"+to_string(n);
      count_jobs.push_back(jname);
    }
  }

  size_t max_cj = _nb_procs;
  FBasedSync _count_sync = FBasedSync(count_jobs, e->SYNCHRO_C, END_TEMP_C, jobs, max_cj, _progress, e, 0);

  uint keep_superk = _keep_tmp ? 1 : 0;
  uint curr = 0;
  for (size_t i=0; i<_bank_paths.size(); i++)
  {
    for (size_t n=0; n<_nb_partitions; n++)
    {
      string command = fmt::format(COUNTER_CMD,
        "",
        e->COUNTER_BIN,                 // counter
        _bank_paths[i],                 // -file
        e->DIR,                         // -run-dir
        _k_size,                        // -kmer-size
        _a_min,                         // -abundance-min
        _max_hash,                      // -max-hash
        _mode,                          // -mode
        1,                              // -nb-cores
        n,                              // -part-id
        _hasher,                        // -hasher
        keep_superk,                      // -keep-tmp
        fmt::format(e->LOG_COUNTER, i, n)  // &> counter.log
      );
      _f_log << command << endl;
      std::system(command.c_str());
      ++*jobs;
      _count_sync.push(_bank_paths[i]+"_"+to_string(n));
      _count_sync.wait_ressources(curr);
      curr++;
    }
  }
  _count_sync.wait_end();
  delete jobs;
  std::system(fmt::format(RM, e->SYNCHRO_C).c_str());
}

void Kmtricks::km_merger()
{
  uint rmemMg = 100;
  size_t max_cj = min(_nb_procs, _max_memory/rmemMg);
  _progress = new ProgressSynchro(
    createIteratorListener(_nb_partitions, "km_merge_within_partition"), System::thread().newSynchronizer()
  );
  
  vector<string> parts;
  for (int p=0; p<_nb_partitions; p++)
    parts.push_back(to_string(p));
  
  size_t* jobs = new size_t(0);
  FBasedSync _merge_sync = FBasedSync(parts, e->SYNCHRO_M, END_TEMP_M, jobs, max_cj, _progress, e, !_keep_tmp);

  string line;
  ifstream fof(_fof_path);
  char cop[256];
  for (size_t i=0; i<_nb_partitions; i++)
  {
    string part_dir = e->STORE_KMERS + fmt::format(PART_DIR, i);
    string path = part_dir + "/partition" + to_string(i) + ".fof";
    ofstream copyfof(path);
    while( getline(fof, line) )
    {
      strcpy(cop, line.c_str());
      string file = basename(cop);
      copyfof << part_dir + "/" + file + ".kmer" << "\n";
    }
    copyfof.close();
    fof.clear();
    fof.seekg(0, ios::beg);

    string command = fmt::format(MERGER_CMD,
      "",
      e->MERGER_BIN,
      e->DIR,         // -run-dir
      i,              // -part-id
      _a_min,         // -abundance-min
      _r_min,         // -recurrence-min
      _mat_str,       // -mode
      fmt::format(e->LOG_MERGER, i)
    );
    _f_log << command << endl;
    std::system(command.c_str());
    ++*jobs;
    _merge_sync.push(parts[i]);
    _merge_sync.wait_ressources(i);
  }
  _merge_sync.wait_end();

  fof.close();
  delete jobs;
  std::system(fmt::format(RM, e->SYNCHRO_M).c_str());
}

void Kmtricks::km_output()
{
  size_t nb_files = _bank_paths.size();
  _progress = new ProgressSynchro(
    createIteratorListener(nb_files, "km_output_convert"), System::thread().newSynchronizer()
  );
  _progress->init();

  string command = fmt::format(OUTPUT_CMD,
    "",
    e->OUTPUT_BIN, // km_output_convert
    e->DIR,        // -run-dir
    nb_files,     // -nb_files
    _str_split,     // -split
    _k_size,       // -kmer-size
    e->LOG_SPLIT
  );
  _f_log << command << endl;
  std::system(command.c_str());
  string odir;
  if (_str_split == "howde") odir = e->STORE_HOWDE;
  else odir = e->STORE_SDSL;

  while (!System::file().doesExist(e->SYNCHRO_SP + END_TEMP_SP))
  {
    sleep(1);
    _progress->set(System::file().listdir(odir).size()-2);
  }
  _progress->finish();
  delete _progress;
}

int main(int argc, char* argv[])
{
  try 
  {
    string mode = "";
    bool env = false;
    if (argc > 1)
      mode = argv[1];
    if (mode == "env")
    {
      env = true; argc--; argv++;
    }
    Kmtricks(env).run(argc, argv);
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
