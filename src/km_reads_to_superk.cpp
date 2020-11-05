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
#include <libgen.h>
#include "km_reads_to_superk.hpp"
#include "signal_handling.hpp"

template<size_t span> struct Functor
{
  void operator() (Parameter parameter)
  {
    KmSuperK&       superk = parameter.superk;
    IProperties*  props = parameter.props;

    IBank* bank = Bank::open(props->getStr(STR_URI_FILE));
    LOCAL (bank);

    char cop[256];
    strcpy(cop, props->getStr(STR_URI_FILE).c_str()); // needed with osx as basename needs a non const char*
    string prefix = basename(cop);
    string _run_dir = props->getStr(STR_RUN_DIR);

    Env* e = new Env(_run_dir, "");
    IteratorListener* _progress(new ProgressSynchro (new IteratorListener(), System::thread().newSynchronizer()));

    _progress->init();
    _progress->setMessage("Compute super-k-mers");

    Storage* _config_storage = StorageFactory(STORAGE_FILE).load(e->STORE_CONFIG);
    Storage* _repart_storage = StorageFactory(STORAGE_FILE).load(e->STORE_PART);

    Configuration _config = Configuration();
    _config.load(_config_storage->getGroup(CONFIG_GRP));
    Repartitor* repartitor = new Repartitor(_repart_storage->getGroup(REPART_GRP));

    string name = e->STORE_SUPERK + fmt::format(TEMP_S, prefix);
    bool lz4 = props->getInt(STR_LZ4_OUT);
    SuperKmerBinFiles* _superKstorage = new SuperKmerBinFiles(name, "superKparts", _config._nb_partitions, lz4);

    typedef typename Kmer<span>::Type                                       Type;
    typedef typename Kmer<span>::ModelCanonical                             ModelCanonical;
    typedef typename Kmer<span>::template ModelMinimizer <ModelCanonical>   Model;

    uint32_t* freq_order = NULL;
    Model model(_config._kmerSize, _config._minim_size, 
      typename gatb::core::kmer::impl::Kmer<span>::ComparatorMinimizerFrequencyOrLex(), freq_order);

    Iterator<Sequence>* itSeq = bank->iterator();
    LOCAL (itSeq);

    BankStats _bankStats;
    PartiInfo<5> pInfo (_config._nb_partitions, _config._minim_size);

    uint _nb_cores = props->getInt(STR_NB_CORES);
    Dispatcher dispatcher(_nb_cores);

    size_t groupSize = 1000;
    bool deleteSync = true;
    dispatcher.iterate(
      itSeq,
      FillPartitions<span, true> (
          model, 1, 0, _config._nb_partitions, _config._nb_cached_items_per_core_per_part,
          _progress, _bankStats, nullptr, *repartitor, pInfo, _superKstorage),
      groupSize, deleteSync
    );

    itSeq->finalize();
    _superKstorage->flushFiles();
    _superKstorage->closeFiles();

    _superKstorage->saveInfoFile(name);
    pInfo.saveInfoFile(name);
    _progress->finish();

    string end_sign = e->SYNCHRO_S + fmt::format(END_TEMP_S, prefix);
    IFile* sync_file = System::file().newFile(end_sign, "w");
    sync_file->flush();
    delete e;
  }
};

KmSuperK::KmSuperK() : Tool ("km_superk")
{
  setParser(new OptionsParser("km_reads_to_superk"));

  getParser()->push_back(new OptionOneParam(STR_URI_FILE,
    "path to read file", true));
  getParser()->push_back(new OptionOneParam(STR_KMER_SIZE,
    "size of a k-mer", true));
  getParser()->push_back(new OptionOneParam(STR_RUN_DIR,
    "kmtricks run directory", true));
  getParser()->push_back(new OptionOneParam(STR_LZ4_OUT,
    "compress output super-k-mers files with lz4 compression", false, "0"));
  getParser()->push_back(new OptionOneParam(STR_NB_CORES,
    "number of cores", true));
}

void KmSuperK::execute()
{
  size_t kmer_size = getInput()->getInt(STR_KMER_SIZE);
  Integer::apply<Functor, Parameter> (kmer_size, Parameter(*this, getInput()));
}

int main(int argc, char* argv[])
{
  try 
  {
    INIT_SIGN;

    KmSuperK().run(argc, argv);
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