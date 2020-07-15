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

#include <fmt/format.h>
#include <tuple>
#include <gatb/kmer/impl/HashSorting.hpp>
#include <gatb/kmer/impl/PartitionsCommand.hpp>
#include <libgen.h>
#include "kmcount.hpp"
#include "CountProcessorDump.hpp"

#ifndef KMTYPE
  typedef uint64_t kmtype_t;
#endif

#ifndef CNTYPE
  typedef uint8_t cntype_t;
#endif

//template <size_t span = KMER_DEFAULT_SPAN>
//
//class CountProcessorDumpPart : public CountProcessorAbstract<span>
//{
//public:
//  typedef typename Kmer<span>::Count Count;
//  typedef typename Kmer<span>::Type Type;
//
//  CountProcessorDumpPart(
//      size_t kmerSize,
//      CountNumber min_abundance,
//      string out_part,
//      uint partId,
//      size_t nbPartsPerPass = 0)
//
//      : _kmerSize(kmerSize), _nbPartsPerPass(nbPartsPerPass), _min_abundance(min_abundance), _out_part(out_part), _partId(partId)
//  {
//    cout << out_part << endl;
//    strcpy(_head, "kmerPart");
//    sprintf(_b, "%03d", _partId);
//    strcat(_head, _b);
//    _part_file.rdbuf()->pubsetbuf(_buffer, 8192);
//    _part_file.open(_out_part, std::ios::app | std::ios::binary);
//    if (!_part_file)
//    {
//      cout << "Unable to open " + _out_part << endl;
//      exit(EXIT_FAILURE);
//    }
//    _part_file.write(_head, strlen(_head));
//
//    model = new typename Kmer<span>::ModelCanonical(kmerSize);
//  }
//
//  virtual ~CountProcessorDumpPart()
//  {
//    _part_file.close();
//  }
//
//  void begin(const Configuration &config)
//  {
//    _nbPartsPerPass = config._nb_partitions;
//    size_t nbTotalPartitions = config._nb_partitions * config._nb_passes;
//  }
//
//  CountProcessorAbstract<span> *clone()
//  {
//    return new CountProcessorDumpPart(_kmerSize, _nbPartsPerPass, _out_part, _partId);
//  }
//
//  void finishClones(vector<ICountProcessor<span> *> &clones)
//  {
//    for (size_t i = 0; i < clones.size(); i++)
//    {
//      if (CountProcessorDumpPart *clone = dynamic_cast<CountProcessorDumpPart *>(clones[i]))
//      {
//        for (map<string, size_t>::iterator it = clone->_namesOccur.begin(); it != clone->_namesOccur.end(); ++it)
//        {
//          this->_namesOccur[it->first] += it->second;
//        }
//      }
//    }
//  }
//
//  void beginPart(size_t passId, size_t partId, size_t cacheSize, const char *name)
//  {
//    size_t actualPartId = partId * (passId * _nbPartsPerPass);
//    _namesOccur[name]++;
//  }
//
//  void endPart(size_t passId, size_t partId)
//  {
//  }
//
//  bool process(size_t partId, const Type &kmer, const CountVector &count, CountNumber sum)
//  {
//    CountNumber kmer_count = count[0];
//    _hk = kmer.getVal();
//    if (kmer_count >= _min_abundance)
//    {
//      _hcount = kmer_count >= 0xFF ? 0xFF : (uint8_t)kmer_count;
//      _part_file.write((char *)&_hk, sizeof(u_int64_t));
//      _part_file.write((char *)&_hcount, sizeof(u_int8_t));
//    }
//    return true;
//  }
//
//private:
//  typename Kmer<span>::ModelCanonical *model;
//  size_t    _kmerSize;
//  size_t    _min_abundance;
//  size_t    _nbPartsPerPass;
//  string    _out_part;
//  ofstream  _part_file;
//  char      _buffer[8192];
//  char      _b[4];
//  char      _head[12];
//  uint64_t  _hk;
//  uint8_t   _hcount;
//  uint      _partId;
//  map<string, size_t> _namesOccur;
//};

template <size_t span>
struct Functor
{
  void operator()(Parameter parameter)
  {
    KmCount &counter = parameter.counter;
    IProperties *props = parameter.props;

    string _run_dir = props->getStr(STR_RUN_DIR);
    uint   hash_mode = props->getInt(STR_MODE);
    uint   nbCores = props->getInt(STR_NB_CORES);
    
    Env *e = new Env(_run_dir, "");
    Storage *_config_storage = StorageFactory(STORAGE_FILE).load(e->STORE_CONFIG);
    Configuration _config = Configuration();
    _config.load(_config_storage->getGroup(CONFIG_GRP));

    uint kmerSize = _config._kmerSize;
    uint max_memory = _config._max_memory;

    char cop[256];
    strcpy(cop, props->getStr(STR_URI_FILE).c_str()); // needed with osx as basename needs a non const char*
    string prefix = basename(cop);
    uint min_abundance = props->getInt(STR_KMER_ABUNDANCE_MIN);

    vector<size_t> nbItemsPerBankPerPart; // vector of length > 1: offsets for multibank counting
    u_int64_t mem = (max_memory * MBYTE) / nbCores;
    typedef typename Kmer<span>::Count Count;
    size_t cacheSize = std::min((u_int64_t)(200 * 1000), mem / (50 * sizeof(Count)));

    string name = e->STORE_SUPERK + fmt::format(TEMP_S, prefix);
    SuperKmerBinFiles *_superKstorage = new SuperKmerBinFiles(name); // a custom constructor that's not in vanilla gatb-core
    PartiInfo<5> pInfo(name);

    uint nb_partitions = _superKstorage->nbFiles();
    e->build_p(nb_partitions - 1);
    MemAllocator pool(nbCores);
    const u_int64_t MBYTE = (1ULL << 20);
    uint64_t memoryPoolSize = max_memory * MBYTE;
    if (pool.getCapacity() == 0)
    {
      pool.reserve(memoryPoolSize);
    }
    else if (memoryPoolSize > pool.getCapacity())
    {
      pool.reserve(0);
      pool.reserve(memoryPoolSize);
    }

    uint part_id = props->getInt(STR_PART_ID);
    uint64_t memReq = ((pInfo.getNbKmer(part_id)*64)/8) + 4096;
    uint64_t hash_mem = 0;
    if (memReq > pool.getCapacity())
    {
      hash_mem = pool.getCapacity();
    }

    uint pass = 0; // probably shouldn't change

    gatb::core::tools::misc::impl::TimeInfo _fillTimeInfo;

    gatb::core::tools::dp::IteratorListener *_progress(new ProgressSynchro(
        new IteratorListener(),
        System::thread().newSynchronizer()));
    _progress->init();

    uint64_t _max_hash = props->getInt(STR_MAX_HASH);
    uint64_t window_size = (uint64_t)ceil((double)_max_hash / (double)nb_partitions);
    
    vector<ICommand*> cmds;
    string path = e->STORE_KMERS + fmt::format(PART_TEMP_K, part_id, prefix);
    CountProcessorDumpPart<span>* dumper = new CountProcessorDumpPart<span>(kmerSize, min_abundance, path, part_id, nb_partitions);

    string hasher = props->getStr(STR_HASHER);
    bool sabuhash = false;
    if (hasher.find("sabuhash") != string::npos) sabuhash = true;

    if (hash_mode)
    {
      if (!hash_mem)
        cmds.push_back(new HashSortingCommand<span>(
          dumper, cacheSize, _progress, _fillTimeInfo, pInfo, pass,
          part_id, 1, kmerSize, pool, nbItemsPerBankPerPart, _superKstorage, window_size, sabuhash));
      else
        cmds.push_back(new HashByHashCommand<span>(
          dumper, cacheSize, _progress, _fillTimeInfo, pInfo, pass,
          part_id, 1, kmerSize, pool, hash_mem, _superKstorage, window_size, sabuhash));
    }
    else
    {
      if (!hash_mem)
        cmds.push_back(new PartitionsByVectorCommand<span>(
          dumper, cacheSize, _progress, _fillTimeInfo, pInfo, pass,
          part_id, 1, kmerSize, pool, nbItemsPerBankPerPart, _superKstorage));
      else
        cmds.push_back(new PartitionsByHashCommand<span>(
          dumper, cacheSize, _progress, _fillTimeInfo, pInfo, pass,
          part_id, 1, kmerSize, pool, hash_mem, _superKstorage));
    }

    cmds[0]->execute();
    string sign = prefix+"_"+to_string(part_id);
    string end_sign = e->SYNCHRO_C + fmt::format(END_TEMP_C, sign);
    IFile *sync_file = System::file().newFile(end_sign, "w");
    sync_file->flush();
    
    System::file().remove(e->STORE_SUPERK + "/" + prefix + ".superk" + "/superKparts." + to_string(part_id));
    pool.free_all();
    _superKstorage->closeFiles();
    _progress->finish();
    delete sync_file;
    delete e;
    delete dumper;
  }
};

KmCount::KmCount() : Tool("km_count")
{
  setParser(new OptionsParser("Kmtricks sub-program: counter"));
  getParser()->setHelp("WARNING: this is a sub-program used by Kmtricks, don't use it directly.");
  getParser()->push_back(new OptionOneParam(STR_URI_FILE,
    "path"));
  getParser()->push_back(new OptionOneParam(STR_RUN_DIR,
    "run directory", true));
  getParser()->push_back(new OptionOneParam(STR_KMER_ABUNDANCE_MIN,
    "abundance min to keep a k-mer", true));
  getParser()->push_back(new OptionOneParam(STR_KMER_SIZE,
    "size of a k-mer", true));
  getParser()->push_back(new OptionOneParam(STR_MAX_HASH,
    "max hash value", false, "0"));
  getParser()->push_back(new OptionOneParam(STR_NB_CORES,
    "nb cores", true));
  getParser()->push_back(new OptionOneParam(STR_MODE,
    "0: k-mers, 1: hashes", false, "1"));
  getParser()->push_back(new OptionOneParam(STR_PART_ID,
    "part id", true));
  getParser()->push_back(new OptionOneParam(STR_HASHER,
    "hash function: [sabuhash || xor]", false, "None"));
}

void KmCount::execute()
{
  /** we get the kmer size chosen by the end user. */
  size_t kmerSize = getInput()->getInt(STR_KMER_SIZE);

  /** We launch dsk with the correct Integer implementation according to the choosen kmer size. */
  Integer::apply<Functor, Parameter>(kmerSize, Parameter(*this, getInput()));
}

/********************************************************************************/

int main(int argc, char *argv[])
{
  try
  {
    KmCount().run(argc, argv);
  }

  catch (OptionFailure &e)
  {
    return e.displayErrors(std::cout);
  }

  catch (Exception &e)
  {
    cerr << "EXCEPTION: " << e.getMessage() << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
