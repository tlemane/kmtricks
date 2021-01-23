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

// TODO: refactoring is needed, ugly duplicate code for testing

#include "km_output_convert.hpp"
#include "signal_handling.hpp"

#include <kmtricks/logging.hpp>
#include <kmtricks/lz4_stream.hpp>
#include <kmtricks/utilities.hpp>
#include <fmt/format.h>

km::log_config km::LOG_CONFIG;

KmConvert::KmConvert(const string &mode) 
  : Tool("km_output_convert"),
    _e(nullptr),
    _howde(false),
    _sdsl(false),
    _vlen(0),
    _filter_size(0),
    _nb_parts(0),
    _win_size(0),
    _kmer_size(0)
{
  if (mode == "from_merge")
  {
    setParser(new OptionsParser("km_output_convert from_merge"));
    getParser()->push_back(new OptionOneParam(STR_RUN_DIR, "kmtricks runtime directory", true));
    getParser()->push_back(new OptionOneParam(STR_NB_FILE, "number of reads files", true));
    getParser()->push_back(new OptionOneParam(STR_SPLIT, "output format: sdsl, howde", true));
    getParser()->push_back(new OptionOneParam(STR_KMER_SIZE, "size of a k-mer", true));
    getParser()->push_back(new OptionOneParam(STR_NB_CORES, "unused but needed by gatb args parser", false, "1"), 0UL, false);

    _from_merge = true;
  }
  else if (mode == "from_count")
  {
    setParser(new OptionsParser("km_output_convert from_count"));
    getParser()->push_back(new OptionOneParam(STR_RUN_DIR, 
      "kmtricks runtime directory", true));
    getParser()->push_back(new OptionOneParam(STR_URI_FILE,
      "file prefix", true));
    getParser()->push_back(new OptionOneParam(STR_KMER_SIZE,
      "size of a k-mer", true));
    getParser()->push_back(new OptionOneParam(STR_SPLIT,
      "output format: sdsl, howde", true));
    getParser()->push_back(new OptionOneParam(STR_NB_CORES,
      "unused but needed by gatb args parser", false, "1"), 0UL, false);

    _from_merge = false;
  }

  km::LOG_CONFIG.show_labels=true;
  km::LOG_CONFIG.level=km::INFO;
}

void KmConvert::parse_args()
{
  _run_dir = getInput()->getStr(STR_RUN_DIR);
  _split_str = getInput()->getStr(STR_SPLIT);
  _howde = filter_format.at(_split_str) != 1;
  _sdsl = true;
  _kmer_size = getInput()->getInt(STR_KMER_SIZE);
  
  if (_from_merge)
  {
    _nb_files = getInput()->getInt(STR_NB_FILE);
    _vlen = NMOD8(NBYTE(_nb_files));
  }
  else
  {
    _f_basename = getInput()->getStr(STR_URI_FILE);
  }
}

void KmConvert::init()
{
  _e = new Env(_run_dir, "");
  if (_from_merge)
  {
    fof_t fof = parse_km_fof(_e->FOF_FILE);
    for (auto& elem: fof)
    {
      string opath;
      if ( _howde )
        opath = _e->STORE_HOWDE + "/" + get<0>(elem) + ".bf";
      else
        opath = _e->STORE_SDSL + "/" + get<0>(elem) + ".sdsl";
      _f_names.push_back(opath);
    }
  }
  ifstream hw(_e->HASHW_MAP, ios::binary | ios::in);

  string mpath;
  uint64_t h0, h1;

  hw.read((char*)&_nb_parts, sizeof(uint32_t));
  for ( int i = 0; i < _nb_parts; i++ )
  {
    hw.read((char *) &h0, sizeof(uint64_t));
    hw.read((char *) &h1, sizeof(uint64_t));
    _hash_windows.push_back(make_tuple(h0, h1));
    if (_from_merge)
    {
      mpath = _e->STORE_MATRIX + fmt::format(BF_T_TEMP, i, i);
      _matrices.push_back(new km::BitMatrixFile<km::IN, km::matrix_t::BF>(mpath));
    }
  }
  _win_size = get<1>(_hash_windows[0]) + 1;
  _filter_size = get<1>(_hash_windows.back()) + 1;

  _sync = _e->SYNCHRO_SP + END_TEMP_SP;
}

void KmConvert::from_merge()
{
  for (int i=0; i<_nb_files; i++)
  {
    bitvector b(_filter_size, 0);
    char* data = (char*)b.data();
    ofstream bf(_f_names[i]);
    uint64_t offset = 0;
    for (int p=0; p<_nb_parts; p++)
    {
      _matrices[p]->read<char>(&data[offset], _win_size/8);
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
      std::memset(header, 0, header_size);
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
      header->numBits      = (uint64_t)_filter_size;
      header->numVectors   = (int)1;
      header->setSizeKnown = false;
      header->setSize      = 0;

      header->info[0].compressor = bvcomp_uncompressed;
      header->info[0].name = 0;
      header->info[0].offset = bw;

      size_t numBytes = sdsl::serialize(b, bf);
      //bw += numBytes;
      header->info[0].numBytes = numBytes;
      header->info[0].filterInfo = (uint64_t)0;

      bf.seekp(ios::beg);
      bf.write((char*)header, header_size);
      bf.close();
      delete[] header;
    }
  }
}

void KmConvert::from_count()
{
  string output_path;
  if (_howde)
    output_path = _e->STORE_HOWDE + "/" + _f_basename + ".bf";
  else
    output_path = _e->STORE_SDSL + "/" + _f_basename + ".sdsl";
  
  ofstream out(output_path, ios::binary | ios::out);
  string in_path;
  string ext = ".kmer.vec";
  bool islz4 = false;
  string path_template = _e->STORE_KMERS + "/partition_{}/{}";
  if (!ifstream(fmt::format(path_template, 0, _f_basename+ext)).good())
  {
    ext = ".kmer.vec.lz4";
    islz4 = true;
  }
  bitvector b(_filter_size, 0);
  char* data = (char*)b.data();
  uint64_t offset = 0;
  for (int i=0; i<_nb_parts; i++)
  {
    in_path = fmt::format(path_template, i, _f_basename+ext);
    km::BitVectorFile<km::IN> in(in_path);
    in.read(&data[offset], _win_size/8);
    offset += _win_size/8;
  }

  if (!_howde)
    sdsl::serialize(b, out);
  else
  {
    uint64_t bneeded = bffileheader_size(1);
    bneeded = round_up_16(bneeded);
    uint32_t header_size = (uint32_t) bneeded;
    bffileheader* header = (bffileheader*)new char[header_size];
    memset(header, 0, header_size);
    header->magic = bffileheaderMagicUn;
    header->headerSize = (uint32_t) sizeof(bffileprefix);
    out.write((char*)header, header_size);
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
    header->numBits      = (uint64_t)_filter_size;
    header->numVectors   = (int)1;
    header->setSizeKnown = false;
    header->setSize      = 0;

    header->info[0].compressor = bvcomp_uncompressed;
    header->info[0].name = 0;
    header->info[0].offset = bw;

    size_t numBytes = sdsl::serialize(b, out);
    //bw += numBytes;
    header->info[0].numBytes = numBytes;
    header->info[0].filterInfo = (uint64_t)0;

    out.seekp(ios::beg);
    out.write((char*)header, header_size);
    out.close();
    delete[] header;
  }
}

void KmConvert::execute()
{
  parse_args();
  init();

  if (_from_merge)
    from_merge();
  else
    from_count();
  
  km::LOG(km::INFO) << "File: " << (_from_merge ? _fof : _f_basename);
  km::LOG(km::INFO) << "Mode: " << (_howde ? "howde" : "sdsl");
  km::LOG(km::INFO) << "Size: " << _filter_size;
  
  for (int i=0; i<_matrices.size(); i++)
    delete _matrices[i];
  delete _e;
}

int main(int argc, char *argv[])
{
  try
  {
    INIT_SIGN;
    
    string mode;
    if (argc > 1)
      mode = argv[1];
    argc--; argv++;

    if (mode == "from_merge" or mode == "from_count")
      KmConvert(mode).run(argc, argv);
    else
    {
      cerr << "km_output_convert subcommands:" << endl;
      cerr << "       km_output_convert from_merge --help" << endl;
      cerr << "       km_output_convert from_count --help" << endl;
      return EXIT_FAILURE;
    }
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
