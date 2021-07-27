/*****************************************************************************
 *   kmtricks
 *   Authors: T. Lemane, R. Chikhi, GATB Team
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
#include <gatb/gatb_core.hpp>
#include <kmtricks/io/kmer_file.hpp>
#include <kmtricks/io/hash_file.hpp>
#include <kmtricks/io/vector_file.hpp>
#include <kmtricks/io/kff_file.hpp>
#include <kmtricks/utils.hpp>
#include <kmtricks/histogram.hpp>

namespace km {

template<size_t span>
class IHashProcessor
{
public:
  virtual bool process(size_t partId, uint64_t hash, const uint32_t count) = 0;
  virtual void finish() = 0;
  virtual ~IHashProcessor() {}
};

template<size_t span>
class ICountProcessor
{
  using Type = typename ::Kmer<span>::Type;
public:
  virtual bool process(size_t partId, const Type& kmer, uint32_t count) = 0;
  virtual void finish() {};
  virtual ~ICountProcessor() {}
};

template<size_t span, size_t MAX_C, size_t buf_size = 32768>
class HashCountProcessor : public IHashProcessor<span>
{
public:
  using Count = typename ::Kmer<span>::Count;
  using Type = typename ::Kmer<span>::Type;
  using km_count_type = typename selectC<DMAX_C>::type;

  HashCountProcessor(uint32_t kmer_size, uint32_t abundance_min, hw_t<MAX_C, buf_size> writer, hist_t hist)
    : m_kmer_size(kmer_size), m_abundance_min(abundance_min), m_writer(writer), m_hist(hist)
  {}

  bool process(size_t partId, uint64_t hash, uint32_t count) override
  {
    if (m_hist) m_hist->inc(count);
    if (count >= m_abundance_min)
    {
      m_count = count >= m_max_c ? m_max_c : static_cast<km_count_type>(count);
      m_writer->write(hash, m_count);
    }
    return true;
  }

  void finish() override { m_writer->flush(); }

private:
  uint32_t m_kmer_size;
  uint32_t m_abundance_min;
  hw_t<MAX_C, buf_size> m_writer;
  hist_t m_hist;
  typename ::Kmer<span>::ModelCanonical m_model {m_kmer_size};
  km_count_type m_count;
  uint32_t m_max_c {std::numeric_limits<km_count_type>::max()};
};

template<size_t span, size_t buf_size = 8192>
class HashVecProcessor : public IHashProcessor<span>
{
public:
  using Count = typename ::Kmer<span>::Count;
  using Type = typename ::Kmer<span>::Type;
  using km_count_type = typename selectC<DMAX_C>::type;

  HashVecProcessor(uint32_t kmer_size, uint32_t abundance_min, bvw_t<buf_size> writer,
                   hist_t hist, size_t window)
    : m_kmer_size(kmer_size), m_abundance_min(abundance_min), m_writer(writer), m_hist(hist),
      m_window(window)
  {
    m_vec.resize(NBYTES(m_window), 0);
  }

  bool process(size_t partId, uint64_t hash, uint32_t count) override
  {
    if (m_hist) m_hist->inc(count);
    if (count >= m_abundance_min)
      BITSET(m_vec, hash - (m_window * partId));
    return true;
  }

  void finish() override { m_writer->write(m_vec); m_writer->flush(); }

private:
  uint32_t m_kmer_size;
  uint32_t m_abundance_min;
  bvw_t<buf_size> m_writer;
  hist_t m_hist;
  typename ::Kmer<span>::ModelCanonical m_model {m_kmer_size};
  km_count_type m_count;
  uint32_t m_max_c {std::numeric_limits<km_count_type>::max()};
  std::vector<uint8_t> m_vec;
  size_t m_window;
};

template<size_t span, size_t MAX_C, size_t buf_size = 8192>
class KmerCountProcessor : public ICountProcessor<span>
{
public:
  using Count = typename ::Kmer<span>::Count;
  using Type = typename ::Kmer<span>::Type;
  using km_count_type = typename selectC<MAX_C>::type;

  KmerCountProcessor(uint32_t kmer_size,
                     uint32_t abundance_min, kw_t<buf_size> writer, hist_t hist)
    : m_kmer_size(kmer_size), m_abundance_min(abundance_min), m_writer(writer), m_hist(hist)
  {}

  bool process(size_t partId, const Type &kmer, uint32_t count) override
  {
    if (m_hist) m_hist->inc(count);
    if (count >= m_abundance_min)
    {
      //Kmer<span> kmkmer(m_model.toString(kmer));
      m_count = count >= m_max_c ? m_max_c : static_cast<km_count_type>(count);
      //m_writer->template write<span, MAX_C>(kmkmer, m_count);
      m_writer->template write_raw<MAX_C>(kmer.get_data(), m_count);
    }
    return true;
  }

private:
  uint32_t m_kmer_size;
  uint32_t m_abundance_min;
  kw_t<buf_size> m_writer;
  hist_t m_hist;
  typename ::Kmer<span>::ModelCanonical m_model {m_kmer_size};
  km_count_type m_count;
  uint32_t m_max_c {std::numeric_limits<km_count_type>::max()};
};

template<size_t span, size_t MAX_C>
class KffCountProcessor : public ICountProcessor<span>
{
public:
  using Count = typename ::Kmer<span>::Count;
  using Type = typename ::Kmer<span>::Type;
  using km_count_type = typename selectC<MAX_C>::type;

  KffCountProcessor(uint32_t kmer_size,
                    uint32_t abundance_min, kff_w_t<MAX_C> writer, hist_t hist)
    : m_kmer_size(kmer_size), m_abundance_min(abundance_min), m_writer(writer), m_hist(hist)
  {}

  bool process(size_t partId, const Type &kmer, uint32_t count) override
  {
    if (m_hist) m_hist->inc(count);
    if (count >= m_abundance_min)
    {
      Kmer<span> kmkmer(m_model.toString(kmer));
      m_count = count >= m_max_c ? m_max_c : static_cast<km_count_type>(count);
      m_writer->template write<span>(kmkmer, m_count);
    }
    return true;
  }

private:
  uint32_t m_kmer_size;
  uint32_t m_abundance_min;
  kff_w_t<MAX_C> m_writer;
  hist_t m_hist;
  typename::Kmer<span>::ModelCanonical m_model {m_kmer_size};
  km_count_type m_count;
  uint32_t m_max_c {std::numeric_limits<km_count_type>::max()};
};

};