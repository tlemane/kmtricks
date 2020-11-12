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

#include <tuple>
#include <libgen.h>
#include <fmt/format.h>
#include <gatb/kmer/impl/HashSorting.hpp>
#include <gatb/kmer/impl/PartitionsCommand.hpp>
#include <kmtricks/logging.hpp>
#include "CountProcessorDump.hpp"
#include "km_superk_to_kmer_counts.hpp"
#include "signal_handling.hpp"

km::log_config km::LOG_CONFIG;

#define NMOD8(byte) ((byte)+(8-((byte)%8)))

template <size_t span>
struct Functor
{
  void operator()(Parameter parameter)
  {
    KmCount &counter = parameter.counter;
    IProperties *props = parameter.props;

    string _run_dir   = props->getStr(STR_RUN_DIR);
    uint   hash_mode  = props->getInt(STR_MODE);
    uint   nbCores    = props->getInt(STR_NB_CORES);
    bool   keep_tmp   = props->getInt(STR_KEEP_TMP);
    bool   lz4        = props->getInt(STR_LZ4_OUT);

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
    uint64_t window_size = NMOD8((uint64_t)ceil((double)_max_hash / (double)nb_partitions));

    vector<ICommand*> cmds;
    string path = e->STORE_KMERS + fmt::format(PART_TEMP_K, part_id, prefix);
    uint64_t vec_size = props->getInt(STR_VEC_ONLY) ? window_size : 0;
    CountProcessorDumpPart<span>* dumper = new CountProcessorDumpPart<span>(kmerSize, min_abundance, path, part_id, lz4, nb_partitions, vec_size);

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

    km::LOG(km::INFO) << "File: " << prefix;
    km::LOG(km::INFO) << "Mode: " << (hash_mode ? "hash" : "kmer");
    km::LOG(km::INFO) << "Out:  " << (vec_size ? "vector" : "value:count");
    km::LOG(km::INFO) << "lz4:  " << (lz4 ? "true" : "false");
    km::LOG(km::INFO) << "Algo: " << (hash_mem ? "ByHash" : "ByVector");

    cmds[0]->execute();
    string sign = prefix+"_"+to_string(part_id);
    string end_sign = e->SYNCHRO_C + fmt::format(END_TEMP_C, sign);
    IFile *sync_file = System::file().newFile(end_sign, "w");
    sync_file->flush();

    if (!keep_tmp)
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
  setParser(new OptionsParser("kmtricks: km_superk_to_kmer_counts"));

  IOptionsParser *hParser = new OptionsParser("hash, only with -mode 1");
  hParser->push_back(new OptionOneParam(STR_HASHER,
    "hash function: sabuhash, xor", false, "xor"));
  hParser->push_back(new OptionOneParam(STR_MAX_HASH,
    "max hash value", false, "0"));

  getParser()->push_back(new OptionOneParam(STR_URI_FILE,
    "path to read file", true));
  getParser()->push_back(new OptionOneParam(STR_RUN_DIR,
    "kmtricks run directory", true));
  getParser()->push_back(new OptionOneParam(STR_KMER_ABUNDANCE_MIN,
    "abundance min to keep a k-mer", false, "2"));
  getParser()->push_back(new OptionOneParam(STR_KMER_SIZE,
    "size of a k-mer", false, "31"));
  getParser()->push_back(new OptionOneParam(STR_PART_ID,
    "partition id", true));
  getParser()->push_back(new OptionOneParam(STR_MODE,
    "0: k-mers, 1: hashes", false, "0"));
  getParser()->push_back(new OptionOneParam(STR_VEC_ONLY,
    "0: hash/count, 1: bit-vector -> when merge is not required", false, "0"));
  getParser()->push_back(new OptionOneParam(STR_NB_CORES,
    "not used, needed by gatb args parser", true));
  getParser()->push_back(new OptionOneParam(STR_KEEP_TMP,
    "keep superkmers files", false, "0"));
  getParser()->push_back(new OptionOneParam(STR_LZ4_OUT,
    "compress output k-mers files with lz4 compression", false, "0"));

  getParser()->push_back(hParser);

  km::LOG_CONFIG.show_labels=true;
  km::LOG_CONFIG.level=km::INFO;
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
    INIT_SIGN;
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
