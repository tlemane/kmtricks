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

#define NMOD8(byte) ((byte)+(8-((byte)%8)))

void wait_end_signal(string sign)
{
  const struct timespec t[] {{0, 500000000L}};
  while (!System::file().doesExist(sign))
  {
    nanosleep(t, nullptr);
  }
}

void signal_callback(int signum)
{
  cout << "\nInterrupt signal (" << signum << ") received" << endl;
  std::system(KILLALL.c_str());
  cout << "All children are killed." << endl;
  exit(EXIT_FAILURE);
}

Kmtricks::Kmtricks()
  : Tool("Kmtricks"),
  _min_hash(0)
{
  getParser()->push_back(new OptionOneParam(STR_URI_FILE,
    "fof that contains path of read files, one per line", true));
  getParser()->push_back(new OptionOneParam(STR_RUN_DIR,
    "directory to write tmp and output files", true));
  getParser()->push_back(new OptionOneParam(STR_KMER_SIZE,
    "size of a kmer", false, "31"));
  getParser()->push_back(new OptionOneParam(STR_KMER_ABUNDANCE_MIN,
    "min abundance threshold for solid kmers", false, "2"));
  getParser()->push_back(new OptionOneParam(STR_KMER_ABUNDANCE_MAX,
    "max abundance threshold for solid kmers", false, "27777777"));
  getParser()->push_back(new OptionOneParam(STR_REC_MIN,
    "min recurrence threshold through datasets for solid kmers", false, "2"));
  getParser()->push_back(new OptionOneParam(STR_NB_PROC,
    "max simultaneously process", false, "8"));
  getParser()->push_back(new OptionOneParam(STR_MAX_MEMORY, 
    "max memory available in megabytes", true));
  getParser()->push_back(new OptionOneParam(STR_UP_TO,
    "run until step N. 1:partitioner, 2:superk, 3:count, 4:merge", false, "0"));
  getParser()->push_back(new OptionOneParam(STR_ONLY,
    "run only step N. 1:partitioner, 2:superk, 3:count, 4:merger", false, "0"));

  IOptionsParser* mParser = new OptionsParser ("kmer count, advanced performance tweaks");

  mParser->push_back (new OptionOneParam (STR_MINIMIZER_TYPE,
    "minimizer type (0=lexi, 1=freq)", false, "0"));
  mParser->push_back (new OptionOneParam (STR_MINIMIZER_SIZE,
    "size of a minimizer", false, "10"));
  mParser->push_back (new OptionOneParam (STR_REPARTITION_TYPE,
    "minimizer repartition (0=unordered, 1=ordered)", false, "0"));
  mParser->push_back(new OptionOneParam ("-nb-parts", "number of partitions", false, "0"));
  getParser()->push_back (mParser);
  
  IOptionsParser* oParser = new OptionsParser("output configuration");

  oParser->push_back(new OptionOneParam(STR_COUNT_SIZE,
    "size of count: 8, 16, 32 bits", false, "8"));
  oParser->push_back(new OptionOneParam(STR_HASHER,
    "hash function: [sabuhash || xor]", false, "xor"));
  oParser->push_back(new OptionOneParam(STR_MAX_HASH,
    "max hash value ( 0 < hash < max(int64) )", false, "1e9"));
  oParser->push_back(new OptionOneParam(STR_MAT_FMT,
    "output matrix format: [0, 1, 2, 3, 4] - 3 and 4 available only with -hasher\n\n"
    "                           0: ascii    <kmer/hash>  <count_0> ... <count_n>\n"
    "                           1: count    <uint64>     <uintX_0> ... <uintX_n>\n"
    "                           2: pa       <uint64>     <presence/absence vecs>\n\n"
    "                           3: bf       (hash0)      <presence/absence vecs>\n"
    "                                          |                    |           \n"
    "                                       (hashN)      <presence/absence vecs>\n\n"
    "                           4: bf_trp   transposition of bf\n", true));
  oParser->push_back(new OptionOneParam(STR_SPLIT,
    "split matrix in indvidual files: [sdsl || howde], (only with -matrix-fmt 4)\n\n"
    "                           sdsl: dump bit vector files with sdsl:bit_vector compatibility\n"
    "                           howde: dump bit vector files with HowDeSBT compatibility", true));
  getParser()->push_back(oParser);
}

void Kmtricks::parse_args()
{
  _fof_path         = getInput()->getStr(STR_URI_FILE);
  _nb_cores         = getInput()->getInt(STR_NB_CORES);
  _max_memory       = getInput()->getInt(STR_MAX_MEMORY);
  _k_size           = getInput()->getInt(STR_KMER_SIZE);
  _a_min            = getInput()->getInt(STR_KMER_ABUNDANCE_MIN);
  _a_max            = getInput()->getInt(STR_KMER_ABUNDANCE_MAX);
  _r_min            = getInput()->getInt(STR_REC_MIN);
  _dir              = getInput()->getStr(STR_RUN_DIR);

  _max_hash         = getInput()->getInt(STR_MAX_HASH);
  _min_hash         = getInput()->getInt(STR_MAX_HASH);
  _mat_fmt          = getInput()->getInt(STR_MAT_FMT);
  _hasher           = getInput()->getStr(STR_HASHER);

  _upto             = getInput()->getInt(STR_UP_TO);
  _only             = getInput()->getInt(STR_ONLY);
  _nb_partitions    = getInput()->getInt("-nb-parts");

  getInput()->add(1, STR_MAX_DISK, "0"); 
  getInput()->add(1, STR_STORAGE_TYPE, "file");
  getInput()->add(1, STR_SOLIDITY_KIND, "sum");
}

void Kmtricks::init()
{
  _nb_procs = _nb_cores;

  char* buffer = new char[1024]; // 

#if __APPLE__
  uint32_t size;
  _NSGetExecutablePath(buffer, &size);
#else
  readlink("/proc/self/exe", buffer, 1024);
#endif

  _path_binary = dirname(buffer);
  delete[] buffer;

  e = new Env(_dir, _path_binary);
  e->build();
  _f_log = ofstream(e->LOG_CMD, ios::out);

  string line, in;
  ifstream fof(_fof_path);
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
  Storage *_config_storage = StorageFactory(STORAGE_FILE).create(e->STORE_CONFIG, true, false);
  ConfigurationAlgorithm<span> configA (bank, getInput());
  configA.execute();
  Configuration _config = configA.getConfiguration();
  if (_nb_partitions != 0)
    _config._nb_partitions = _nb_partitions;
  _config.save(_config_storage->getGroup("config"));
  _nb_partitions = _config._nb_partitions;
  e->build_p(_nb_partitions);

  auto window_size = NMOD8((uint64_t)ceil((double)_max_hash / (double)_nb_partitions));

  for (int i=0; i<_nb_partitions; i++)
  {
    _hash_windows.emplace_back(i*window_size, ((i+1)*window_size)-1);
  }  
  
  size_t maxCores   = _nb_cores;
  size_t maxMem     = _max_memory;
  size_t minMemJob  = 1000;
  
  size_t maxJ_core  = min(maxCores/2, _bank_paths.size());
  maxJ_core         = max(maxJ_core, (size_t)1);
  size_t maxJ_mem   = maxMem/minMemJob;
  maxJ_mem          = max(maxJ_mem, (size_t)1);
  size_t maxJobs    = min(maxJ_core, maxJ_mem);
  _max_job          = maxJobs;

  _max_job          = max(_max_job, maxCores);
  _core_p_j         = maxCores/_max_job;
  _core_p_j         = max(_core_p_j, (size_t)1);

  _mem_p_j          = maxMem / _max_job;
  _mem_p_j          = max(_mem_p_j, (size_t)minMemJob);
  

}

void Kmtricks::execute()
{
  parse_args();
  init();
  signal(SIGINT, signal_callback);

  bool all = (_only == 0);

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

  end:
    _f_log.close();
    cout << endl;
}

void Kmtricks::km_part()
{
  string command = fmt::format(PARTITIONER_CMD,
    "",
    e->PARTITIONER_BIN, // partitioner
    _fof_path,          // -file
    _k_size,            // -kmer-size
    _nb_procs,          // -nb-procs
    e->DIR,             // -run-dir
    e->SYNCHRO_P,       // -dir-synchro
    e->LOG_PARTITIONER  // &> partitioner.log
  );
  _f_log << command << endl;
  std::system(command.c_str());
  wait_end_signal(e->SYNCHRO_P + END_TEMP_P);
}

void Kmtricks::km_superk()
{

  uint rmemSk = 300;
  size_t max_cj = min(_nb_procs, (_max_memory/rmemSk));

  _progress = new ProgressSynchro(
    createIteratorListener (_bank_paths.size()*_nb_partitions, "km_superk"), System::thread().newSynchronizer());

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
      e->SYNCHRO_S,                 // -dir-synchro
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
}

void Kmtricks::km_count()
{
  uint64_t nbs = _bank_paths.size()*_nb_partitions;
  _progress = new ProgressSynchro(
    createIteratorListener (nbs, "km_count"), System::thread().newSynchronizer());

  size_t* jobs = new size_t(0);
  vector<string> count_jobs;
  //vector<uint64_t> vmem(nbs);
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
      uint64_t nbelem = tmpinfo.getNbKmer(n);
      //vmem.push_back((nbelem*64)/8000000);
      jname = _bank_paths[i]+"_"+to_string(n);
      count_jobs.push_back(jname);
    }
  }

  //uint rmemCt = *max_element(vmem.begin(), vmem.end()) + 50;
  size_t max_cj = _nb_procs;//min(_nb_procs, (_max_memory/rmemCt));
  FBasedSync _count_sync = FBasedSync(count_jobs, e->SYNCHRO_C, END_TEMP_C, jobs, max_cj, _progress, e, 0);
  
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
}

void Kmtricks::km_merger()
{
  uint rmemMg = 100;
  size_t max_cj = min(_nb_procs, _max_memory/rmemMg);
  _progress = new ProgressSynchro(
    createIteratorListener(_nb_partitions, "km_merge"), System::thread().newSynchronizer()
  );
  
  vector<string> parts;
  for (int p=0; p<_nb_partitions; p++)
    parts.push_back(to_string(p));
  
  size_t* jobs = new size_t(0);
  FBasedSync _merge_sync = FBasedSync(parts, e->SYNCHRO_M, END_TEMP_M, jobs, max_cj, _progress, e, 0);
  for (size_t i=0; i<_nb_partitions; i++)
  {
    string part_dir = e->STORE_KMERS + fmt::format(PART_DIR, i);
    string fof = part_dir + "/partition" + to_string(i) + ".fof";
    vector<string> partition_paths = System::file().listdir(part_dir);
    ofstream out_fof(fof, ios::out);
    for (auto& path: partition_paths)
    {
      if ( (path.size()>2) && path.find("fof") == string::npos)
      {
        out_fof << e->STORE_KMERS + fmt::format(PART_DIR, i) + "/" + path + "\n";
      }
    }

    string command = fmt::format(MERGER_CMD,
      "",
      e->MERGER_BIN,  // kbf_merger
      fof,            // -file
      e->DIR,         // -run-dir
      get<0>(_hash_windows[i]),      // -min-hash
      get<1>(_hash_windows[i]),      // -max-hash
      i,              // -part-id
      _a_min,         // -abundance-min
      _r_min,         // -recurrence-min
      _mat_fmt,       // -mode
      fmt::format(e->LOG_MERGER, i)
    );
    _f_log << command << endl;
    std::system(command.c_str());
    ++*jobs;
    _merge_sync.push(parts[i]);
    _merge_sync.wait_ressources(i);
  }
  _merge_sync.wait_end();

  delete jobs;
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
