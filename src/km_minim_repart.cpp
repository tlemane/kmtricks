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

#include <utility>
#include <stdexcept>
#include <exception>
#include <kmtricks/logging.hpp>
#include <kmtricks/utilities.hpp>
#include "km_minim_repart.hpp"
#include "signal_handling.hpp"

km::log_config km::LOG_CONFIG;

Repart::Repart() : Tool ("km_part")
{
  setParser(new OptionsParser("km_minim_repart"));

  getParser()->push_back(new OptionOneParam(STR_URI_FILE,
                                              "fof that contains path of read files, one per line", true));
  getParser()->push_back(new OptionOneParam(STR_KMER_SIZE,
                                              "size of a k-mer", true));
  getParser()->push_back(new OptionOneParam(STR_RUN_DIR,
                                              "kmtricks run directory", true));
   getParser()->push_back(new OptionOneParam(STR_NB_CORES,
                                              "number of cores", false, "8"));

  km::LOG_CONFIG.show_labels=true;
  km::LOG_CONFIG.level=km::INFO;
}

template<size_t span> struct Functor 
{ 
  void operator() (Parameter parameter)
  {
    Repart&      rep   = parameter.rep;
    IProperties* props = parameter.props;
    
    Env* e = new Env(props->getStr(STR_RUN_DIR), "");
    fof_t fof = parse_km_fof(props->getStr(STR_URI_FILE));
    string in = all_files(fof);
    km::LOG(km::INFO) << in;
    IBank *bank = Bank::open(in);
    LOCAL (bank);
 
    Storage* _repart_storage = StorageFactory(STORAGE_FILE).create(
        e->STORE_PART,
        true,
        false);
    
    Storage* _config_storage = StorageFactory(STORAGE_FILE).load(
        e->STORE_CONFIG);
    
    Configuration _config = Configuration();
    _config.load(_config_storage->getGroup(CONFIG_GRP));

    RepartitorAlgorithm<span> repart(
        bank,
        _repart_storage->getGroup(REPART_GRP),
        _config,
        props->getInt(STR_NB_CORES));
    repart.execute ();

    km::LOG(km::INFO) << "Repartition file at " << e->STORE_PART;
    
    IFile* sync_file =  System::file().newFile(e->SYNCHRO_P+END_TEMP_P, "w");
    sync_file->flush();
    delete sync_file;
    delete e;
  }
};

void Repart::execute()
{
  size_t kmer_size = getInput()->getInt(STR_KMER_SIZE);
  Integer::apply<Functor, Parameter> (kmer_size, Parameter(*this, getInput()));
}

int main(int argc, char *argv[])
{
  try 
  {
    INIT_SIGN;

    Repart().run(argc, argv);
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