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

#include <kmtricks/logging.hpp>
#include "km_merge_within_partition.hpp"
#include "signal_handling.hpp"

km::log_config km::LOG_CONFIG;

KmMerge::KmMerge() : Tool("km_merge")
{
  setParser(new OptionsParser("km_merge_within_partition"));

  getParser()->push_back(new OptionOneParam(STR_RUN_DIR, "kmtricks run directory", true));
  getParser()->push_back(new OptionOneParam(STR_PART_ID, "partition id", true));
  getParser()->push_back(new OptionOneParam(STR_KMER_ABUNDANCE_MIN, "abundance min to keep a k-mer", true));
  getParser()->push_back(new OptionOneParam(STR_REC_MIN, "recurrence min to keep a k-mer", true));
  getParser()->push_back(new OptionOneParam(STR_SAVE_IF, "save a non-solid k-mer if it occurs in N other datasets", false, "0"));
  getParser()->push_back(new OptionOneParam(STR_MODE, "output matrix format: ascii, bin, pa, bf, bf_trp"));
  getParser()->push_back(new OptionOneParam(STR_HSIZE, "file header size in byte", false, "0"));
  getParser()->push_back(new OptionOneParam(STR_NB_CORES, "not used, needed by gatb args parser", false, "1"));

  km::LOG_CONFIG.show_labels=true;
  km::LOG_CONFIG.level=km::INFO;
}

void KmMerge::parse_args()
{
  _run_dir        = getInput()->getStr(STR_RUN_DIR);
  _min_r          = getInput()->getInt(STR_REC_MIN);
  _id             = getInput()->getInt(STR_PART_ID);
  _mode           = output_format.at(getInput()->getStr(STR_MODE));
  
  string tmp      = getInput()->getStr(STR_KMER_ABUNDANCE_MIN);
  if (System::file().doesExist(tmp))
  {
    ifstream in(tmp, ios::in);
    string line;
    while (getline(in, line))
      _abs_vec.push_back(stoi(line));
    _min_a = 0;
  }
  else
    _min_a          = getInput()->getInt(STR_KMER_ABUNDANCE_MIN);
  
  _save_if        = getInput()->getInt(STR_SAVE_IF);

  e = new Env(_run_dir, "");
  _fofpath = fmt::format(e->STORE_KMERS + PART_DIR + "/partition{}.fof", _id, _id);

  ifstream hw(e->HASHW_MAP, ios::in);
  uint64_t h0, h1;
  uint32_t nb_parts;
  hw.read((char*)&nb_parts, sizeof(uint32_t));
  for ( int i = 0; i < nb_parts; i++ )
  {
    hw.read((char *) &h0, sizeof(uint64_t));
    hw.read((char *) &h1, sizeof(uint64_t));
    _hash_windows.push_back(make_tuple(h0, h1));
  }
  hw.close();
  _lower_hash = get<0>(_hash_windows[_id]);
  _upper_hash = get<1>(_hash_windows[_id]);
}

void KmMerge::merge_to_pa_matrix()
{
  ofstream fout;
  string opath = e->STORE_MATRIX + fmt::format(PA_TEMP, _id, _id);
  fout.open(opath, ios::binary);
  while (!_m->end)
  {
    _m->next();
    if (_m->keep)
    {
      fout.write((char*)&_m->m_khash, sizeof(kmtype_t));
      fout.write((char*)_m->_bit_vector, _m->vlen);
    }
  }
  fout.close();
  string end_sign = e->SYNCHRO_M + fmt::format(END_TEMP_M, _id);
  IFile *sync_file = System::file().newFile(end_sign, "w");
  sync_file->flush();
  delete sync_file;
}

void KmMerge::merge_to_bin()
{
  ofstream fout;
  string opath = e->STORE_MATRIX + fmt::format(CO_TEMP, _id, _id);
  fout.open(opath, ios::binary | ios::out);
  while (!_m->end)
  {
    _m->next();
    if (_m->keep)
    {
      fout.write((char*)&_m->m_khash, sizeof(kmtype_t));
      for (size_t i=0; i<_m->nb_files; i++)
      {
        fout.write((char*)&_m->counts[i], sizeof(cntype_t));
      }
    }
  }
  fout.close();
  string end_sign = e->SYNCHRO_M + fmt::format(END_TEMP_M, _id);
  IFile *sync_file = System::file().newFile(end_sign, "w");
  sync_file->flush();
  delete sync_file;
}

void KmMerge::merge_to_ascii()
{
  string opath = e->STORE_MATRIX + fmt::format(AS_TEMP, _id, _id);
  ofstream fout;
  fout.open(opath, ios::out);
  while (!_m->end)
  {
    _m->next();
    if (_m->keep)
    {
#if KTYPE == 128
      fout << _m->m_khash;
#else
      fout << to_string(_m->m_khash);
#endif
      for (size_t vi=0; vi<_m->nb_files; vi++)
        fout << " " << to_string(_m->counts[vi]);
      fout << "\n";
    }
  }
  fout.close();
  string end_sign = e->SYNCHRO_M + fmt::format(END_TEMP_M, _id);
  IFile *sync_file = System::file().newFile(end_sign, "w");
  sync_file->flush();
  delete sync_file;
}

void KmMerge::merge_to_bf_pa()
{
  ofstream fout;
  string opath = e->STORE_MATRIX + fmt::format(BF_NT_TEMP, _id, _id);
  fout.open(opath, ios::binary);
  uchar* empty = new uchar[_m->vlen]();
  uint64_t current_hash = _lower_hash;

  while(!_m->end)
  {
    _m->next();
    while(_m->m_khash > current_hash)
    {
      fout.write((char*)empty, _m->vlen);
      current_hash++;
    }

    if (_m->keep)
    {
      fout.write((char*)_m->_bit_vector, _m->vlen);
      current_hash = _m->m_khash+1;
    }
  }

  while (current_hash <= _upper_hash)
  {
    fout.write((char*)empty, _m->vlen);
    current_hash++;
  }
  fout.close();

  string end_sign = e->SYNCHRO_M + fmt::format(END_TEMP_M, _id);
  IFile *sync_file = System::file().newFile(end_sign, "w");
  sync_file->flush();
  delete sync_file;
}

void KmMerge::transpose()
{
  string path_mat = e->STORE_MATRIX + fmt::format(BF_NT_TEMP, _id, _id);
  uint n = _upper_hash - _lower_hash + 1;
  uint m = NMOD8(NBYTE(_m->nb_files));
  BitMatrix *mat = new BitMatrix(path_mat, n, m, true);
  BitMatrix *trp = mat->transpose();
  string outp = e->STORE_MATRIX + fmt::format(BF_T_TEMP, _id, _id);
  remove(path_mat.c_str());
  trp->dump(outp);
  delete mat;
  delete trp;
}

void KmMerge::execute()
{

  parse_args();
  size_t hsize = getInput()->getInt(STR_HSIZE);
  bool setbv = _mode > 1;

  km::LOG(km::INFO) << "Fof:   " << _fofpath;
  km::LOG(km::INFO) << "Mode:  " << output_format_str.at(_mode);
  km::LOG(km::INFO) << "A-min: " << _min_a;
  km::LOG(km::INFO) << "R-min: " << _min_r;
  km::LOG(km::INFO) << "Save-if: " << _save_if;

  if (!_min_a && _abs_vec.size() > 0)
    _m = new Merger<kmtype_t, cntype_t>(
      _fofpath, _abs_vec, _min_r, hsize, setbv, _save_if, true);
  else
    _m = new Merger<kmtype_t, cntype_t>(
      _fofpath, _min_a, _min_r, hsize, setbv, _save_if, true);
  
  switch (_mode)
{
  case 0:
    merge_to_ascii();
    break;
  case 1:
    merge_to_bin();
    break;
  case 2:
    merge_to_pa_matrix();
    break;
  case 3:
    merge_to_bf_pa();
    break;
  case 4:
    merge_to_bf_pa();
    transpose();

  default:
    break;
  }

  km::LOG(km::INFO) << "ABS VEC: " << _abs_vec;
  km::LOG(km::INFO) << "NON_SOLID: " << _m->_non_solid;
  km::LOG(km::INFO) << "SAVED: " << _m->_saved;
  delete _m;
  delete e;
}

int main (int argc, char* argv[])
{
  try
  {
    INIT_SIGN;

    KmMerge().run(argc, argv);
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
