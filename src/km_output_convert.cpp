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
#include "km_output_convert.hpp"

KmConvert::KmConvert() : Tool("km_output_convert")
{
  setParser(new OptionsParser("km_output_convert"));
  getParser()->push_back(new OptionOneParam(STR_RUN_DIR, "kmtricks runtime directory", true));
  getParser()->push_back(new OptionOneParam(STR_URI_FILE, "extract bitvectors only for the files contained in the fof", false, "0"));
  getParser()->push_back(new OptionOneParam(STR_NB_FILE, "number of reads files", true));
  getParser()->push_back(new OptionOneParam(STR_NB_PARTS, "number of partitions", true));
  getParser()->push_back(new OptionOneParam(STR_SPLIT, "output format: sdsl, howde", true));
  getParser()->push_back(new OptionOneParam(STR_KMER_SIZE, "size of a k-mer", true));
  getParser()->push_back(new OptionOneParam(STR_VERBOSE, "verbosity level", false, "1"));
  getParser()->push_back(new OptionOneParam(STR_NB_CORES, "unused but needed by gatb args parser", false, "1"));
}

void KmConvert::parse_args()
{
  _run_dir = getInput()->getStr(STR_RUN_DIR);
  _split_str = getInput()->getStr(STR_SPLIT);
  _sdsl = true;
  _howde = filter_format.at(_split_str) != 1;
  _nb_files = getInput()->getInt(STR_NB_FILE);
  _vlen = NMOD8(NBYTE(_nb_files));
  _nb_parts = getInput()->getInt(STR_NB_PARTS);
  _kmer_size = getInput()->getInt(STR_KMER_SIZE);
  _fof = getInput()->getStr(STR_URI_FILE);
}

void KmConvert::init()
{
  _e = new Env(_run_dir, "");
  ifstream fof(_e->FOF_FILE);
  string line;
  char cop[256];
  string fname;
  string opath;
  uint cnt = 0;
  while ( getline(fof, line))
  {
    _fof_pos[line] = cnt;
    cnt++;
    strcpy(cop, line.c_str());
    fname = basename(cop);
    if ( _howde )
      opath = _e->STORE_HOWDE + "/" + fname + ".bf";
    else
      opath = _e->STORE_SDSL + "/" + fname + ".sdsl";
    _f_names.push_back(opath);
  }
  fof.close();

  ifstream hw(_e->HASHW_MAP, ios::binary | ios::in);

  string mpath;
  uint64_t h0, h1;
  for ( int i = 0; i < _nb_parts; i++ )
  {
    hw.read((char *) &h0, sizeof(uint64_t));
    hw.read((char *) &h1, sizeof(uint64_t));
    _hash_windows.push_back(make_tuple(h0, h1));
    mpath = _e->STORE_MATRIX + fmt::format(BF_T_TEMP, i, i);
    _matrices.push_back(ifstream(mpath, ios::in | ios::binary));
  }
  _win_size = get<1>(_hash_windows[0]) + 1;
  _filter_size = get<1>(_hash_windows.back()) + 1;

  if ( _fof != "0" )
  {
    ifstream infof(_fof);
    string f;
    while (getline(infof, f))
      _pos.push_back(_fof_pos[f]);
    infof.close();
  }
  _sync = _e->SYNCHRO_SP + END_TEMP_SP;
  delete _e;
}
void KmConvert::execute()
{
  parse_args();
  init();
  for (int i=0; i<_nb_files; i++)
  {
    if ( _fof != "0")
      if (!(std::find(_pos.begin(), _pos.end(), i) != _pos.end()))
        continue;

    bitvector b(_filter_size, 0);
    char* data = (char*)b.data();
    ofstream bf(_f_names[i]);
    uint64_t offset = 0;
    for (int p=0; p<_nb_parts; p++)
    {
      _matrices[p].read((char*)&data[offset], _win_size/8);
      offset += _win_size/8;
    }
    if (!_howde)
      sdsl::serialize(b, bf);
    else
    {
      uint64_t bneeded = bffileheader_size(1);
      bneeded = round_up_16(bneeded);
      uint32_t header_size = (uint32_t) bneeded;
      bffileheader* header = (bffileheader*)new char[header_size];
      memset(header, 0, header_size);
      header->magic = bffileheaderMagicUn;
      header->headerSize = (uint32_t) sizeof(bffileprefix);
      bf.write((char*)header, header_size);
      size_t bw = header_size;

      header->magic        = bffileheaderMagic;
      header->headerSize   = header_size;
      header->version      = bffileheaderVersion;
      header->bfKind       = bfkind_simple;
      header->padding1     = 0;
      header->kmerSize     = _kmer_size;
      header->numHashes    = (uint32_t)1;
      header->hashSeed1    = (uint64_t)0;
      header->hashSeed2    = (uint64_t)0;
      header->hashModulus  = (uint64_t)_filter_size;
      header->numBits      = (uint64_t)2;
      header->numVectors   = (int)1;
      header->setSizeKnown = false;
      header->setSize      = 0;

      header->info[0].compressor = bvcomp_uncompressed;
      header->info[0].name = 0;
      header->info[0].offset = bw;

      size_t numBytes = sdsl::serialize(b, bf);
      bw += numBytes;
      header->info[0].numBytes = numBytes;
      header->info[0].filterInfo = (uint64_t)0;

      bf.seekp(ios::beg);
      bf.write((char*)header, header_size);
      bf.close();
      delete[] header;
    }
  }
  IFile* sync_file =  System::file().newFile(_sync, "w");
  sync_file->flush();
  delete sync_file;
}

int main(int argc, char *argv[])
{
  try
  {
    KmConvert().run(argc, argv);
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
