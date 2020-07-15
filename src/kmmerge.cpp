#include "kmmerge.hpp"

KmMerge::KmMerge() : Tool("km_merge")
{
  setParser(new OptionsParser("Kmtricks sub-program: merger"));

  getParser()->push_back(new OptionOneParam(STR_URI_FILE, "fof file", true));
  getParser()->push_back(new OptionOneParam(STR_RUN_DIR, "run directory", true));
  getParser()->push_back(new OptionOneParam(STR_MIN_HASH, "lower bound hash", true));
  getParser()->push_back(new OptionOneParam(STR_MAX_HASH, "upper bound hash", true));
  getParser()->push_back(new OptionOneParam(STR_PART_ID, "partition id", true));
  getParser()->push_back(new OptionOneParam(STR_KMER_ABUNDANCE_MIN, "abundance min", true));
  getParser()->push_back(new OptionOneParam(STR_REC_MIN, "recurrence min", true));
  getParser()->push_back(new OptionOneParam(STR_MODE, "output format: [0, 1, 2, 3, 4], see kmtricks help"));
  getParser()->push_back(new OptionOneParam(STR_HSIZE, "header size in byte", false, "12"));
  getParser()->push_back(new OptionOneParam(STR_NB_CORES, "nb cores", false, "0"));
}

void KmMerge::parse_args()
{
  _run_dir        = getInput()->getStr(STR_RUN_DIR);
  _fofpath        = getInput()->getStr(STR_URI_FILE);
  _min_a          = getInput()->getInt(STR_KMER_ABUNDANCE_MIN);
  _min_r          = getInput()->getInt(STR_REC_MIN);
  _lower_hash     = getInput()->getInt(STR_MIN_HASH);
  _upper_hash     = getInput()->getInt(STR_MAX_HASH);
  _id             = getInput()->getInt(STR_PART_ID);
  _mode           = getInput()->getInt(STR_MODE);
  e = new Env(_run_dir, "");
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
      fout.write((char*)&_m->m_khash, sizeof(KType));
      fout.write((char*)_m->_bit_vector, _m->vlen);
    }
  }
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
      fout.write((char*)&_m->m_khash, sizeof(KType));
      for (size_t i=0; i<_m->nb_files; i++)
      {
        fout.write((char*)&_m->counts[i], sizeof(CType));
      }
    }
  }
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
      fout << to_string(_m->m_khash);
      for (size_t vi=0; vi<_m->nb_files; vi++)
        fout << " " << to_string(_m->counts[vi]);
      fout << "\n";
    }
  }
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
  uchar* empty = new uchar[_m->vlen];
  uint64_t current_hash = _lower_hash;
  while (!_m->end)
  {
    _m->next();
    while (_m->m_khash > current_hash)
    {
      fout.write((char*)empty, _m->vlen);
      current_hash++;
    }

    if (_m->keep)
    {
      fout.write((char*)_m->_bit_vector, _m->vlen);
      current_hash = _m->m_khash;
    }
  }
  while (current_hash < _upper_hash)
  {
    fout.write((char*)empty, _m->vlen);
    current_hash++;
  }

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
  trp->dump(outp);
  delete mat;
  delete trp;
}
void KmMerge::execute()
{
  parse_args();
  size_t hsize = getInput()->getInt(STR_HSIZE);
  bool setbv = _mode > 1;

  _m = new Merger<KType, CType>(_fofpath, _min_a, _min_r, hsize, setbv);
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

  _m->destroy();
  delete _m;
}

int main (int argc, char* argv[])
{
  try
  {
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
