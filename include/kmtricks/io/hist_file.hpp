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
#include <kmtricks/io/io_common.hpp>
#include <kmtricks/utils.hpp>
#include <kmtricks/histogram.hpp>

namespace km {

class HistFileHeader : public KmHeader
{
public:
  HistFileHeader() {};

  void serialize(std::ostream* stream) override
  {
    _serialize(stream);
    stream->write(reinterpret_cast<char*>(&hist_magic), sizeof(hist_magic));
    stream->write(reinterpret_cast<char*>(&kmer_size), sizeof(kmer_size));
    stream->write(reinterpret_cast<char*>(&id), sizeof(id));
    stream->write(reinterpret_cast<char*>(&lower), sizeof(lower));
    stream->write(reinterpret_cast<char*>(&upper), sizeof(upper));
    stream->write(reinterpret_cast<char*>(&uniq), sizeof(uniq));
    stream->write(reinterpret_cast<char*>(&total), sizeof(total));
    stream->write(reinterpret_cast<char*>(&oob_ln), sizeof(oob_ln));
    stream->write(reinterpret_cast<char*>(&oob_lu), sizeof(oob_lu));
    stream->write(reinterpret_cast<char*>(&oob_un), sizeof(oob_un));
    stream->write(reinterpret_cast<char*>(&oob_uu), sizeof(oob_uu));
  }

  void deserialize(std::istream* stream) override
  {
    _deserialize(stream);
    stream->read(reinterpret_cast<char*>(&hist_magic), sizeof(hist_magic));
    stream->read(reinterpret_cast<char*>(&kmer_size), sizeof(kmer_size));
    stream->read(reinterpret_cast<char*>(&id), sizeof(id));
    stream->read(reinterpret_cast<char*>(&lower), sizeof(lower));
    stream->read(reinterpret_cast<char*>(&upper), sizeof(upper));
    stream->read(reinterpret_cast<char*>(&uniq), sizeof(uniq));
    stream->read(reinterpret_cast<char*>(&total), sizeof(total));
    stream->read(reinterpret_cast<char*>(&oob_ln), sizeof(oob_ln));
    stream->read(reinterpret_cast<char*>(&oob_lu), sizeof(oob_lu));
    stream->read(reinterpret_cast<char*>(&oob_un), sizeof(oob_un));
    stream->read(reinterpret_cast<char*>(&oob_uu), sizeof(oob_uu));
  }

  void sanity_check() override
  {
    _sanity_check();
    if (hist_magic != MAGICS.at(KM_FILE::HIST))
      throw IOError("Invalid file format.");
  }

public:
  uint64_t hist_magic {MAGICS.at(KM_FILE::HIST)};
  uint32_t kmer_size;
  uint32_t id;
  uint64_t lower;
  uint64_t upper;
  uint64_t uniq;
  uint64_t total;
  uint64_t oob_lu;
  uint64_t oob_uu;
  uint64_t oob_ln;
  uint64_t oob_un;
};

template<size_t buf_size = 8192>
class HistWriter : public IFile<HistFileHeader, std::ostream, buf_size>
{
  using ocstream = lz4_stream::basic_ostream<buf_size>;
public:
  HistWriter(const std::string& path, KHist& hist, bool lz4)
    : IFile<HistFileHeader, std::ostream, buf_size>(path, std::ios::out | std::ios::binary)
  {
    this->m_header.compressed = lz4;
    this->m_header.kmer_size = hist.m_ksize;
    this->m_header.id = hist.m_idx;
    this->m_header.lower = hist.m_lower;
    this->m_header.upper = hist.m_upper;
    this->m_header.uniq = hist.m_uniq;
    this->m_header.total = hist.m_total;
    this->m_header.oob_ln = hist.m_oob_ln;
    this->m_header.oob_lu = hist.m_oob_lu;
    this->m_header.oob_un = hist.m_oob_un;
    this->m_header.oob_uu = hist.m_oob_uu;

    this->m_header.serialize(this->m_first_layer.get());

    this->template set_second_layer<ocstream>(this->m_header.compressed);

    this->m_second_layer->write(reinterpret_cast<char*>(hist.m_hist_u.data()),
                                  hist.m_hist_u.size() * sizeof(uint64_t));
    this->m_second_layer->write(reinterpret_cast<char*>(hist.m_hist_n.data()),
                                  hist.m_hist_n.size() * sizeof(uint64_t));
  }
};

template<size_t buf_size = 8192>
class HistReader : public IFile<HistFileHeader, std::istream, buf_size>
{
  using icstream = lz4_stream::basic_istream<buf_size>;
public:
  HistReader(const std::string& path)
    : IFile<HistFileHeader, std::istream, buf_size>(path, std::ios::in | std::ios::binary)
  {
    this->m_header.deserialize(this->m_first_layer.get());
    this->m_header.sanity_check();

    this->template set_second_layer<icstream>(this->m_header.compressed);
  }

  hist_t get()
  {
    hist_t histo = std::make_shared<KHist>(this->m_header.id, this->m_header.kmer_size,
                                           this->m_header.lower, this->m_header.upper);
    histo->m_oob_lu = this->m_header.oob_lu;
    histo->m_oob_uu = this->m_header.oob_uu;
    histo->m_oob_ln = this->m_header.oob_ln;
    histo->m_oob_un = this->m_header.oob_un;
    this->m_second_layer->read(
      reinterpret_cast<char*>(histo->m_hist_u.data()),
      histo->m_hist_u.size()*sizeof(uint64_t));
    this->m_second_layer->read(
      reinterpret_cast<char*>(histo->m_hist_n.data()),
      histo->m_hist_n.size()*sizeof(uint64_t));
    return histo;
  }

  void write_as_text(std::ostream& stream, bool n)
  {
    hist_t histo = get();
    uint64_t current = histo->lower();
    stream << "@LOWER=" << histo->lower() << "\n";
    stream << "@UPPER=" << histo->upper() << "\n";

    if (n)
    {
      histo->set_type(KHistType::TOTAL);
      stream << "@OOB_L=" << histo->oob_lower_total() << "\n";
      stream << "@OOB_U=" << histo->oob_upper_total() << "\n";
      for_each(histo->begin(), histo->end(), [&current, &stream](uint64_t c) {
        stream << std::to_string(current) << " " << std::to_string(c) << "\n";
        current++;
      });
    }
    else
    {
      stream << "@OOB_L=" << histo->oob_lower_unique() << "\n";
      stream << "@OOB_U=" << histo->oob_upper_unique() << "\n";
      for_each(histo->begin(), histo->end(), [&current, &stream](uint64_t c) {
        stream << std::to_string(current) << " " << std::to_string(c) << "\n";
        current++;
      });
    }
  }
};

};
