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

#include "kmpart.hpp"

string get_str_fof(string fof_path)
{
  string s;
  string line;
  ifstream in_fof(fof_path);

  while ( getline (in_fof, line) ) { s += line + ","; }
  s.pop_back();
  return s;
}

Repart::Repart() : Tool ("km_part")
{
  setParser(new OptionsParser("Kmtricks sub-program: partitioner"));
  getParser()->push_back(new OptionOneParam(STR_URI_FILE,
    "fof that contains one fastx per line", true));
  getParser()->push_back(new OptionOneParam(STR_KMER_SIZE,
    "size of a k-mer", true));
  getParser()->push_back(new OptionOneParam(STR_RUN_DIR,
    "root of run directory", true));
  getParser()->push_back(new OptionOneParam(STR_DIR_SYNCHRO,
    "directory to write synchronization files", true));
  getParser()->push_back(new OptionOneParam(STR_NB_CORES,
    "nb cores"));
};

template<size_t span> struct Functor 
{ 
  void operator() (Parameter parameter)
  {
    Repart&      rep   = parameter.rep;
    IProperties* props = parameter.props;
    
    Env* e = new Env(props->getStr(STR_RUN_DIR), "");
    string in = get_str_fof(props->getStr(STR_URI_FILE));
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