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

#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <map>

#include <gatb/gatb_core.hpp>
#include <kmtricks/lz4_stream.hpp>
#include <kmtricks/utilities.hpp>
#include <kmtricks/io.hpp>
#include <kff_io.hpp>
#include "config.hpp"


template <size_t span = KMER_DEFAULT_SPAN>
class CountProcessorDumpPart : public CountProcessorAbstract<span>
{
public:
  typedef typename Kmer<span>::Count  Count;
  typedef typename Kmer<span>::Type   Type;

  CountProcessorDumpPart(
    size_t kmerSize,
    CountNumber min_abundance,
    const string& out_part,
    uint partId,
    bool lz4,
    km::KmerFile<km::OUT, kmtype_t, cntype_t>* cmf,
    km::KHist *khist,
    size_t nbPartsPerPass = 0,
    size_t window = 0)

    : _kmerSize(kmerSize), _nbPartsPerPass(nbPartsPerPass), _lz4_output(lz4),
      _min_abundance(min_abundance), _out_part(out_part), _partId(partId),
      _window(window/8), _window_bits(window), _hk(0), _hcount(0), _cmf(cmf), _khist(khist)
  {
    if (_window)
    {
      _out_part += ".vec";
      _vec.resize(_window);
    }
    model = new typename Kmer<span>::ModelCanonical(kmerSize);
    _mc = maxc.at(sizeof(cntype_t));

    if (_window)
      _vec.resize(_window);
  }

  CountProcessorDumpPart(
    size_t kmerSize,
    CountNumber min_abundance,
    const string& out_part,
    uint partId,
    bool lz4,
    km::BitVectorFile<km::OUT>* bvf,
    km::KHist* khist,
    size_t nbPartsPerPass = 0,
    size_t window = 0)

    : _kmerSize(kmerSize), _nbPartsPerPass(nbPartsPerPass), _lz4_output(lz4),
      _min_abundance(min_abundance), _out_part(out_part), _partId(partId),
      _window(window/8), _window_bits(window), _hk(0), _hcount(0), _bvf(bvf), _khist(khist)
  {
    if (_window)
    {
      _vec.resize(_window);
    }
    
    model = new typename Kmer<span>::ModelCanonical(kmerSize);
    _mc = maxc.at(sizeof(cntype_t));

    if (_window)
      _vec.resize(_window);
  }

  CountProcessorDumpPart(
    size_t kmerSize,
    CountNumber min_abundance,
    uint partId,
    Kff_file* kff_file,
    km::KHist* khist,
    size_t nbPartsPerPass = 0,
    size_t window = 0
  ) : _kmerSize(kmerSize), _min_abundance(min_abundance), _kff_file(kff_file), _khist(khist),
      _nbPartsPerPass(nbPartsPerPass), _window(window/8), _window_bits(window)
  {
    uint8_t encoding[] = {0, 1, 3, 2};
    _kff_file->write_encoding(encoding);

    Section_GV sgv(_kff_file);
    sgv.write_var("k", _kmerSize);
    sgv.write_var("max", 1);
    sgv.write_var("data_size", sizeof(cntype_t));
    sgv.close();

    _sr = new Section_Raw(_kff_file);

    model = new typename Kmer<span>::ModelCanonical(kmerSize);
  }

  ~CountProcessorDumpPart()
  {
    if (_window)
      flush();
    delete _sr;
    delete model;
  }

  void begin(const Configuration &config)
  {
    _nbPartsPerPass = config._nb_partitions;
    size_t nbTotalPartitions = config._nb_partitions * config._nb_passes;
  }

  CountProcessorAbstract<span> *clone()
  {
    if (_window)
    {
      return new CountProcessorDumpPart(
        _kmerSize, _min_abundance, _out_part, _partId, _lz4_output, _bvf, _khist, _nbPartsPerPass, _window_bits);
    }
    else if (_kff_file)
    {
      return new CountProcessorDumpPart(
        _kmerSize, _min_abundance, _partId, _kff_file, _khist, _nbPartsPerPass, _window_bits
      );
    }
    else
    {
      return new CountProcessorDumpPart(
        _kmerSize, _min_abundance, _out_part, _partId, _lz4_output, _cmf, _khist, _nbPartsPerPass, _window_bits);
    }
  }

  void finishClones(vector<ICountProcessor<span> *> &clones)
  {
    for (size_t i = 0; i < clones.size(); i++)
    {
      if (CountProcessorDumpPart *clone = dynamic_cast<CountProcessorDumpPart *>(clones[i]))
      {
        for (map<string, size_t>::iterator it = clone->_namesOccur.begin(); it != clone->_namesOccur.end(); ++it)
        {
          this->_namesOccur[it->first] += it->second;
        }
      }
    }
  }

  void beginPart(size_t passId, size_t partId, size_t cacheSize, const char *name)
  {
    size_t actualPartId = partId * (passId * _nbPartsPerPass);
    _namesOccur[name]++;
  }

  void endPart(size_t passId, size_t partId)
  { 
    if (_sr)
      _sr->close();
  }

  uint8_t uint8_packing(std::string sequence)
  {
    size_t size = sequence.length();
    assert(size <= 4);

    uint8_t val = 0;
    for (size_t i=0; i<size; i++)
    {
      val <<= 2;
      val += (sequence[i] >> 1) & 0b11;
    }
    return val;
  }

  void encode_sequence(std::string sequence, uint8_t* encoded)
  {
    size_t size = sequence.length();
    size_t remnant = size % 4;
    if (remnant > 0)
    {
      encoded[0] = uint8_packing(sequence.substr(0, remnant));
      encoded += 1;
    }

    size_t nb_uint_needed = size/4;
    for (size_t i=0; i<nb_uint_needed; i++)
    {
      encoded[i] = uint8_packing(sequence.substr(remnant + 4*i, 4));
    }
  }

  void u8from32(uint8_t b[4], uint32_t u32)
  {
    b[3] = (uint8_t)u32;
    b[2] = (uint8_t)(u32>>=8);
    b[1] = (uint8_t)(u32>>=8);
    b[0] = (uint8_t)(u32>>=8);
  }

  void u8from16(uint8_t b[2], uint16_t u16)
  {
    b[1] = (uint8_t)(u16>>=8);
    b[0] = (uint8_t)(u16>>=8);
  }

  bool process(size_t partId, const Type &kmer, const CountVector &count, CountNumber sum)
  {
    CountNumber kmer_count = count[0];
    _khist->inc(kmer_count);
    _hk = kmer.getVal();
    if (kmer_count >= _min_abundance)
    {
      if (_kff_file)
      {
        std::string seq = model->toString(kmer);
        std::cerr << seq << std::endl;
        uint8_t encoded[1024] = {0};
        uint8_t counts[sizeof(cntype_t)];
        
        encode_sequence(seq.c_str(), encoded);
        _hcount = kmer_count >= _mc ? _mc : (cntype_t)kmer_count;
        
        if (sizeof(cntype_t) == 1) counts[0] = _hcount;
        else if (sizeof(cntype_t) == 2) u8from16(counts, _hcount);
        else u8from32(counts, _hcount);
        
        _sr->write_compacted_sequence(encoded, _kmerSize, counts);
      }
      else if (_window)
      {
        BITSET(_vec, _hk -(_window_bits*_partId));
      }
      else
      {
        _hcount = kmer_count >= _mc ? _mc : (cntype_t)kmer_count;
        _cmf->write(_hk, _hcount);
      }
    }
    return true;
  }

  void flush()
  {
    _bvf->write(_vec);
  }

private:
  typename Kmer<span>::ModelCanonical *model;
  size_t    _kmerSize;
  size_t    _min_abundance;
  size_t    _nbPartsPerPass;
  string    _out_part;
  kmtype_t  _hk;
  cntype_t  _hcount;
  uint      _partId;
  bool      _lz4_output;
  size_t    _window;
  size_t    _window_bits;
  map<string, size_t> _namesOccur;
  vector<char>     _vec;
  uint64_t  _mc;
  km::BitVectorFile<km::OUT>* _bvf {nullptr};
  km::KmerFile<km::OUT, kmtype_t, cntype_t>* _cmf {nullptr};
  km::KHist* _khist {nullptr};
  Kff_file* _kff_file {nullptr};
  Section_Raw* _sr {nullptr};
};
