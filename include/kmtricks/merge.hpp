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
#include <vector>
#include <kmtricks/kmer.hpp>
#include <kmtricks/utils.hpp>
#include <kmtricks/io/matrix_file.hpp>
#include <kmtricks/io/pa_matrix_file.hpp>
#include <kmtricks/io/kmer_file.hpp>
#include <kmtricks/io/hash_file.hpp>
#include <kmtricks/io/vector_matrix_file.hpp>
#include <kmtricks/packc.hpp>

#ifdef WITH_PLUGIN
#include <kmtricks/plugin_manager.hpp>
#include <kmtricks/plugin.hpp>
#endif

namespace km {

template<size_t MAX_K, size_t MAX_C>
class IMergeObserver
{
  using count_type = typename selectC<MAX_C>::type;
public:
  IMergeObserver() {}
  virtual void process(Kmer<MAX_K>& kmer, std::vector<count_type>& counts) = 0;
};

template<size_t MAX_K, size_t MAX_C>
using imo_t = std::shared_ptr<IMergeObserver<MAX_K, MAX_C>>;

template<size_t MAX_C>
class MergeStatistics
{
  using count_type = typename selectC<MAX_C>::type;
public:
  MergeStatistics(size_t nb_files)
   : m_nb_files(nb_files)
  {
    m_non_solid.resize(m_nb_files, 0);
    m_rescued.resize(m_nb_files, 0);
    m_uniq_wo_rescue.resize(m_nb_files, 0);
    m_uniq_w_rescue.resize(m_nb_files, 0);
    m_total_wo_rescue.resize(m_nb_files, 0);
    m_total_w_rescue.resize(m_nb_files, 0);
  }

  void inc_ns(uint32_t i) { m_non_solid[i]++; }
  void inc_rd(uint32_t i) { m_rescued[i]++; }
  void inc_uwo(uint32_t i) { m_uniq_w_rescue[i]++; m_uniq_wo_rescue[i]++; }
  void inc_uw(uint32_t i) { m_uniq_w_rescue[i]++; }
  void inc_two(uint32_t i, count_type c) { m_total_wo_rescue[i] += c; m_total_w_rescue[i] += c; }
  void inc_tw(uint32_t i, count_type c) { m_total_w_rescue[i] += c; }

  void serialize(const std::string& path)
  {
    std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
    {
      out << "NON_SOLID" << '\t'; write_vector(out, m_non_solid, '\t') << "\n";
      out << "RESCUED" << '\t'; write_vector(out, m_rescued, '\t') << "\n";
      out << "UNIQUE_WO_RESCUE" << '\t'; write_vector(out, m_uniq_wo_rescue, '\t') << "\n";
      out << "UNIQUE_W_RESCUE" << '\t'; write_vector(out, m_uniq_w_rescue, '\t') << "\n";
      out << "TOTAL_WO_RESCUE" << '\t'; write_vector(out, m_total_wo_rescue, '\t') << "\n";
      out << "TOTAL_W_RESCUE" << '\t'; write_vector(out, m_total_w_rescue, '\t') << "\n";
    }
  }

  const std::vector<uint64_t>& get_non_solid() const { return m_non_solid; }
  const std::vector<uint64_t>& get_rescued() const { return m_rescued; }
  const std::vector<uint64_t>& get_unique_wo_rescue() const { return m_uniq_wo_rescue; }
  const std::vector<uint64_t>& get_unique_w_rescue() const { return m_uniq_w_rescue; }
  const std::vector<uint64_t>& get_total_wo_rescue() const { return m_total_wo_rescue; }
  const std::vector<uint64_t>& get_total_w_rescue() const { return m_total_w_rescue; }

private:
  size_t m_nb_files {0};
  std::vector<uint64_t> m_non_solid;
  std::vector<uint64_t> m_rescued;
  std::vector<uint64_t> m_uniq_wo_rescue;
  std::vector<uint64_t> m_uniq_w_rescue;
  std::vector<uint64_t> m_total_wo_rescue;
  std::vector<uint64_t> m_total_w_rescue;
};

template<size_t MAX_K, size_t MAX_C>
class KmerMerger
{
  using count_type = typename selectC<MAX_C>::type;
public:
  struct element
  {
    Kmer<MAX_K> value;
    count_type count {0};
    bool is_set {false};
  };

public:
  KmerMerger(std::vector<std::string>& paths,
         std::vector<uint32_t>& abundance_min_vec,
         uint32_t kmer_size,
         uint32_t recurrence_min,
         uint32_t save_if)
    : m_paths(paths), m_a_min_vec(abundance_min_vec), m_kmer_size(kmer_size),
      m_r_min(recurrence_min), m_save_if(save_if)
  {
    init_stream();
    init_state();
  }

  const Kmer<MAX_K>& current() const
  {
    return m_current;
  }

  const std::vector<count_type>& counts() const
  {
    return m_counts;
  }

  MergeStatistics<MAX_C>* get_infos() const
  {
    return m_infos.get();
  }

  void init_stream()
  {
    for (auto& path: m_paths)
      m_input_streams.push_back(std::make_shared<KmerReader<8192>>(path));
    m_size = m_paths.size();
    m_kmer_size = m_input_streams[0]->infos().kmer_size;
  }

  void init_state()
  {
    for (size_t i=0; i<m_size; i++)
    {
      m_elements.push_back(element{});
      m_elements[i].value.set_k(m_kmer_size);

      if (read_next(i))
        m_elements[i].is_set = true;

      if ( (!m_current_set || m_elements[i].value < m_current) && m_elements[i].is_set )
      {
        m_current = m_next = m_elements[i].value;
        m_current_set = true;
      }
    }
    m_counts.resize(m_size, 0);
    m_infos = std::make_unique<MergeStatistics<MAX_C>>(m_size);
  }


#ifdef WITH_PLUGIN
  void set_plugin(IMergePlugin* plugin)
  {
    m_plugin = plugin;
  }
#endif

  bool keep()
  {
    return m_keep;
  }

  bool next()
  {
    m_keep = false;
    m_finish = true;
    m_next_set = false;

    uint32_t recurrence = 0;
    uint32_t solid_in = 0;
    m_current = m_next;
    m_need_check.clear();
    for (size_t i=0; i<m_size; i++)
    {
      if (m_elements[i].is_set && m_elements[i].value == m_current)
      {
        m_finish = false;
        m_counts[i] = m_elements[i].count;
        if (m_counts[i] >= m_a_min_vec[i])
        {
          recurrence++;
          solid_in++;

          if (m_infos)
          {
            m_infos->inc_two(i, m_counts[i]);
            m_infos->inc_uwo(i);
          }
        }
        else
        {
          if (m_infos)
            m_infos->inc_ns(i);
          if (m_save_if)
            m_need_check.push_back(i);
          else
            m_counts[i] = 0;
        }
        if (!read_next(i))
          m_elements[i].is_set = false;
      }
      else
      {
        m_counts[i] = 0;
      }

      if (m_elements[i].is_set && (!m_next_set || m_elements[i].value < m_next))
      {
        m_next = m_elements[i].value;
        m_next_set = true;
      }
    }

    for (auto& f : m_need_check)
    {
      if (!(solid_in >= m_save_if))
        m_counts[f] = 0;
      else
      {
        if (m_infos)
        {
          m_infos->inc_rd(f);
          m_infos->inc_uw(f);
          m_infos->inc_tw(f, m_counts[f]);
        }
      }
    }

    if (recurrence >= m_r_min)
      m_keep = true;

#ifdef WITH_PLUGIN
    if (m_plugin)
    {
      m_keep = m_plugin->process_kmer(m_current.get_data64(), m_counts);
    }
#endif

    return !m_finish;
  }

  void write_as_bin(const std::string& path, bool compressed)
  {
    MatrixWriter mw(path, m_kmer_size, 1, m_size, 0, m_partition, compressed);
    while (next())
    {
      if (m_keep)
      {
        mw.template write<MAX_K, MAX_C>(m_current, m_counts);
      }
    }
  }

  void write_as_pa(const std::string& path, bool compressed)
  {
    PAMatrixWriter pw(path, m_kmer_size, m_size, 0, m_partition, compressed);
    std::vector<uint8_t> bit_vec(NBYTES(m_size));
    while (next())
    {
      if (m_keep)
      {
        set_bit_vector(bit_vec, m_counts);
        pw.template write<MAX_K>(m_current, bit_vec);
      }
    }
  }

  void write_as_pa_text(const std::string& path)
  {
    std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
    while (next())
    {
      if (m_keep)
      {
        out << m_current.to_string();
        for (auto& c: m_counts)
          out << " " << (c > 0 ? '1' : '0');
        out << "\n";
      }
    }
  }

  void write_as_text(const std::string& path)
  {
    std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
    while (next())
    {
      if (m_keep)
      {
        out << m_current.to_string();
        for (auto& c: m_counts)
          out << " " << std::to_string(c);
        out << "\n";
      }
    }
  }

  void merge(imo_t<MAX_K, MAX_C> obs)
  {
    while (next())
      if (m_keep)
        obs->process(m_current, m_counts);
  }

private:
  bool read_next(size_t i)
  {
    return m_input_streams[i]->template read<MAX_K, MAX_C>(m_elements[i].value,
                                                           m_elements[i].count);
  }

private:
  std::vector<std::string>& m_paths;
  uint32_t m_a_min;
  uint32_t m_r_min;
  uint32_t m_save_if;
  uint32_t m_partition;

  std::vector<kr_t<8192>> m_input_streams;
  std::vector<element> m_elements;
  std::vector<size_t> m_need_check;

  uint32_t m_size;
  uint32_t m_kmer_size;
  std::vector<uint32_t>& m_a_min_vec;

  Kmer<MAX_K> m_next;
  Kmer<MAX_K> m_current;
  bool m_next_set {false};
  bool m_current_set {false};
  std::vector<count_type> m_counts;

  bool m_keep {false};
  bool m_finish {false};

  std::unique_ptr<MergeStatistics<MAX_C>> m_infos {nullptr};

#ifdef WITH_PLUGIN
  IMergePlugin* m_plugin {nullptr};
#endif
};

template<size_t MAX_C, size_t buf_size = 8192, typename Reader = HashReader<buf_size>>
class HashMerger
{
  using count_type = typename selectC<MAX_C>::type;
public:
  struct element
  {
    uint64_t value;
    count_type count {0};
    bool is_set {false};
  };

public:
  HashMerger(std::vector<std::string>& paths,
         std::vector<uint32_t>& abundance_min_vec,
         uint32_t recurrence_min,
         uint32_t save_if)
    : m_paths(paths), m_a_min_vec(abundance_min_vec),
      m_r_min(recurrence_min), m_save_if(save_if)
  {
    init_stream();
    init_state();
  }

  uint64_t current() const
  {
    return m_current;
  }

  const std::vector<count_type>& counts() const
  {
    return m_counts;
  }

  MergeStatistics<MAX_C>* get_infos() const
  {
    return m_infos.get();
  }

  void init_stream()
  {
    for (auto& path: m_paths)
      m_input_streams.push_back(std::make_shared<Reader>(path));
    m_size = m_paths.size();
    m_partition = m_input_streams[0]->infos().partition;
  }

  void init_state()
  {
    for (size_t i=0; i<m_size; i++)
    {
      m_elements.push_back(element{});

      if (read_next(i))
        m_elements[i].is_set = true;

      if ( (!m_current_set || m_elements[i].value < m_current) && m_elements[i].is_set )
      {
        m_current = m_next = m_elements[i].value;
        m_current_set = true;
      }
    }
    m_counts.resize(m_size, 0);
    m_infos = std::make_unique<MergeStatistics<MAX_C>>(m_size);
  }

#ifdef WITH_PLUGIN
  void set_plugin(IMergePlugin* plugin)
  {
    m_plugin = plugin;
  }
#endif

  bool keep()
  {
    return m_keep;
  }

  bool next()
  {
    m_keep = false;
    m_finish = true;
    m_next_set = false;

    uint32_t recurrence = 0;
    uint32_t solid_in = 0;
    m_need_check.clear();
    m_current = m_next;
    for (size_t i=0; i<m_size; i++)
    {
      if (m_elements[i].is_set && m_elements[i].value == m_current)
      {
        m_finish = false;
        m_counts[i] = m_elements[i].count;

        if (m_counts[i] >= m_a_min_vec[i])
        {
          recurrence++;
          solid_in++;

          if (m_infos)
          {
            m_infos->inc_two(i, m_counts[i]);
            m_infos->inc_uwo(i);
          }
        }
        else
        {
          if (m_infos)
            m_infos->inc_ns(i);
          if (m_save_if)
            m_need_check.push_back(i);
          else
            m_counts[i] = 0;
        }
        if (!read_next(i))
          m_elements[i].is_set = false;
      }
      else
      {
        m_counts[i] = 0;
      }

      if (m_elements[i].is_set && (!m_next_set || m_elements[i].value < m_next))
      {
        m_next = m_elements[i].value;
        m_next_set = true;
      }
    }
    for (auto& f : m_need_check)
    {
      if (!(solid_in >= m_save_if))
        m_counts[f] = 0;
      else
      {
        if (m_infos)
        {
          m_infos->inc_rd(f);
          m_infos->inc_uw(f);
          m_infos->inc_tw(f, m_counts[f]);
        }
      }
    }
    if (recurrence >= m_r_min)
      m_keep = true;

#ifdef WITH_PLUGIN
    if (m_plugin)
    {
      m_keep = m_plugin->process_hash(m_current, m_counts);
    }
#endif

    return !m_finish;
  }

  void write_as_bin(const std::string& path, bool compressed)
  {
    MatrixHashWriter<8192> mhw(path, sizeof(m_counts[0]), m_size, 0, m_partition, compressed);
    while (next())
    {
      if (m_keep)
      {
        mhw.template write<MAX_C>(m_current, m_counts);
      }
    }
  }

  void write_as_text(const std::string& path)
  {
    std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
    while (next())
    {
      if (m_keep)
      {
        out << std::to_string(m_current);
        for (auto& c: m_counts)
          out << " " << std::to_string(c);
        out << "\n";
      }
    }
  }

  void write_as_pa(const std::string& path, bool compressed)
  {
    PAHashMatrixWriter<8192> phw(path, m_size, 0, m_partition, compressed);
    std::vector<uint8_t> bit_vec(NBYTES(m_size));
    while (next())
    {
      if (m_keep)
      {
        set_bit_vector(bit_vec, m_counts);
        phw.write(m_current, bit_vec);
      }
    }
  }

  void write_as_pa_text(const std::string& path)
  {
    std::ofstream out(path, std::ios::out); check_fstream_good(path, out);
    while (next())
    {
      if (m_keep)
      {
        out << std::to_string(m_current);
        for (auto& c: m_counts)
          out << " " << (c > 0 ? '1' : '0');
        out << "\n";
      }
    }
  }

  void write_as_bf(const std::string& path, uint64_t lower, uint64_t upper, bool compressed)
  {
    std::vector<uint8_t> bit_vec(NBYTES(m_size), 0);
    std::vector<uint8_t> empty_vec(NBYTES(m_size), 0);
    uint64_t current = lower;
    VectorMatrixWriter<8192> vmw(path, m_size, 0, m_partition, lower, upper-lower+1, compressed);
    while (next())
    {
      while (m_current > current)
      {
        vmw.write(empty_vec);
        current++;
      }
      if (m_keep)
      {
        set_bit_vector(bit_vec, m_counts);
        vmw.write(bit_vec);
        current = m_current + 1;
      }
    }
    while (current <= upper)
    {
      vmw.write(empty_vec);
      current++;
    }
  }

  void write_as_bfc(const std::string& path, uint64_t lower, uint64_t upper, int w, bool compressed)
  {
    std::vector<uint8_t> cbit_vec(byte_count_pack(m_size, w), 0);
    std::vector<uint8_t> empty_vec(byte_count_pack(m_size, w), 0);
    uint64_t current = lower;

    VectorMatrixWriter<8192> vmw(path, m_size * w, 0, m_partition, lower, upper-lower+1, compressed);

    while (next())
    {
      while (m_current > current)
      {
        vmw.write(empty_vec);
        current++;
      }
      if (m_keep)
      {
        pack_v(m_counts, cbit_vec, w);
        vmw.write(cbit_vec);
        current = m_current + 1;
      }
    }
    while (current <= upper)
    {
      vmw.write(empty_vec);
      current++;
    }
  }

  void write_as_bft(const std::string& path, uint64_t lower, uint64_t upper, bool compressed)
  {
    write_as_bf(path+".tmp", lower, upper, compressed);
    BitMatrix mat(ROUND_UP(upper-lower+1, 8), ROUND_UP(m_size, 8)/8, true);
    {
      VectorMatrixReader vmr(path+".tmp");
      vmr.load(mat);
    }
    std::remove(std::string(path+".tmp").c_str());
    BitMatrix *trp = mat.transpose();
    VectorMatrixWriter<8192> vmw(path, m_size, 0, m_partition, lower, upper-lower+1, compressed);
    vmw.dump(*trp);
    delete trp;
  }

private:
  bool read_next(size_t i)
  {
    if constexpr(std::is_same_v<Reader, HashReader<buf_size>>)
      return m_input_streams[i]->template read<MAX_C>(m_elements[i].value,
                                                      m_elements[i].count);
    else
      return m_input_streams[i]->read(m_elements[i].value,
                                      m_elements[i].count);
  }

private:
  std::vector<std::string>& m_paths;
  uint32_t m_a_min;
  uint32_t m_r_min;
  uint32_t m_save_if;
  uint32_t m_partition;

  std::vector<std::shared_ptr<Reader>> m_input_streams;
  std::vector<element> m_elements;
  std::vector<size_t> m_need_check;

  uint32_t m_size;
  std::vector<uint32_t>& m_a_min_vec;

  uint64_t m_next;
  uint64_t m_current;
  bool m_next_set {false};
  bool m_current_set {false};
  std::vector<count_type> m_counts;

  bool m_keep {false};
  bool m_finish {false};

  std::unique_ptr<MergeStatistics<MAX_C>> m_infos {nullptr};

#ifdef WITH_PLUGIN
  IMergePlugin* m_plugin {nullptr};
#endif
};
};
