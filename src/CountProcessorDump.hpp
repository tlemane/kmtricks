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
    string out_part,
    uint partId,
    bool lz4,
    size_t nbPartsPerPass = 0,
    size_t window = 0)

    : _kmerSize(kmerSize), _nbPartsPerPass(nbPartsPerPass), _lz4_output(lz4), _min_abundance(min_abundance), _out_part(out_part), _partId(partId), _window(window/8), _window_bits(window)
  {
    strcpy(_head, "kmerPart");
    sprintf(_b, "%03d", _partId);
    strcat(_head, _b);
    _part_file.rdbuf()->pubsetbuf(_buffer, 8192);
    
    if (_window) _out_part += ".vec";
    
    if (_lz4_output)
    {
      _part_file.open(_out_part+".lz4", std::ios::app | std::ios::binary);
      _lzstream = new lz4_stream::ostream(_part_file);
      _writer = _lzstream;
    }
    else
    {
      _part_file.open(_out_part, std::ios::app | std::ios::binary);
      _writer = &_part_file;
    }

    if (!_part_file)
    {
      cout << "Unable to open " + _out_part << endl;
      exit(EXIT_FAILURE);
    }
    
    model = new typename Kmer<span>::ModelCanonical(kmerSize);
    _mc = maxc.at(sizeof(cntype_t));

    if (_window)
      _vec.resize(_window);
    else
      _writer->write(_head, strlen(_head));
  }

  ~CountProcessorDumpPart()
  {
    if (_window)
      flush();
    if (_lzstream)
      _lzstream->close();
    _part_file.close();

  }

  void begin(const Configuration &config)
  {
    _nbPartsPerPass = config._nb_partitions;
    size_t nbTotalPartitions = config._nb_partitions * config._nb_passes;
  }

  CountProcessorAbstract<span> *clone()
  {
    return new CountProcessorDumpPart(_kmerSize, _min_abundance, _out_part, _partId, _window, _lz4_output);
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
  }

  bool process(size_t partId, const Type &kmer, const CountVector &count, CountNumber sum)
  {
    CountNumber kmer_count = count[0];
    _hk = kmer.getVal();
    if (kmer_count >= _min_abundance)
    {
      if (_window)
      {
        BITSET(_vec, _hk -(_window_bits*_partId));
      }
      else
      {
        _hcount = kmer_count >= _mc ? _mc : (cntype_t)kmer_count;
        _writer->write((char*)&_hk, sizeof(kmtype_t));
        _writer->write((char *)&_hcount, sizeof(cntype_t));
      }
    }
    return true;
  }

  void flush()
  {
    std::copy(_vec.begin(), _vec.end(), std::ostream_iterator<uint8_t>(*_writer));
  }

private:
  typename Kmer<span>::ModelCanonical *model;
  size_t    _kmerSize;
  size_t    _min_abundance;
  size_t    _nbPartsPerPass;
  string    _out_part;
  ofstream  _part_file;
  ostream   *_writer;
  char      _buffer[8192];
  char      _b[4];
  char      _head[12];
  kmtype_t  _hk;
  cntype_t  _hcount;
  uint      _partId;
  bool      _lz4_output;
  size_t    _window;
  size_t    _window_bits;
  lz4_stream::ostream *_lzstream;
  map<string, size_t> _namesOccur;
  vector<uint8_t>     _vec;
  uint64_t  _mc;
};
