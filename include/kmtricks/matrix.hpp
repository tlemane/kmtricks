#pragma once

#include <vector>
#include <string>
#include <filesystem>
#include <tuple>
#include <memory>

#include <kmtricks/utils.hpp>
#include <kmtricks/io/pa_matrix_file.hpp>
#include <kmtricks/io/kmer_file.hpp>
#include <kmtricks/itask.hpp>
#include <kmtricks/task_pool.hpp>

namespace fs = std::filesystem;

namespace km {

template<std::size_t MAX_K, std::size_t MAX_C>
class FilterTask : public ITask
{
  using count_type = typename selectC<MAX_C>::type;

  public:
    FilterTask(const std::string& matrix,
               const std::string& kmers,
               const std::string& output,
               const std::string& koutput,
               const std::string& vec,
               bool cpr, bool count,
               const std::tuple<bool, bool, bool>& out_types)
      : ITask(0),
        m_matrix(matrix),
        m_kmers(kmers),
        m_output(output),
        m_koutput(koutput),
        m_vec(vec),
        m_cpr(cpr),
        m_count(count),
        m_out_types(out_types)
    {
    }

    void exec()
    {
      if (m_count)
        f_count_matrix(m_out_types);
      else
        f_pa_matrix(m_out_types);
    }

    void preprocess() {}
    void postprocess()
    {
      fs::remove(m_kmers);
    }

  private:

    void f_count_matrix(const std::tuple<bool, bool, bool>& out_types)
    {
      const bool with_vector = std::get<0>(out_types);
      const bool with_matrix = std::get<1>(out_types);
      const bool with_kmer = std::get<2>(out_types);

      std::unique_ptr<std::ofstream> vout {nullptr};
      std::unique_ptr<MatrixWriter<8192>> mw {nullptr};
      std::unique_ptr<KmerWriter<8192>> kw {nullptr};

      KmerReader<8192> kr(m_kmers);
      Kmer<MAX_K> kmer; kmer.set_k(kr.infos().kmer_size);
      count_type count;

      MatrixReader<8192> mr(m_matrix);
      Kmer<MAX_K> kmer2; kmer2.set_k(mr.infos().kmer_size);

      if (with_vector)
      {
        vout = std::make_unique<std::ofstream>(m_vec, std::ios::out);
      }

      if (with_matrix)
      {
        mw = std::make_unique<MatrixWriter<8192>>(m_output,
                                                  mr.infos().kmer_size,
                                                  mr.infos().count_slots,
                                                  mr.infos().nb_counts + 1,
                                                  mr.infos().id,
                                                  mr.infos().partition,
                                                  m_cpr);
      }

      if (with_kmer)
      {
        kw = std::make_unique<KmerWriter<8192>>(m_koutput,
                                                kr.infos().kmer_size,
                                                kr.infos().count_slots,
                                                kr.infos().id,
                                                kr.infos().partition,
                                                m_cpr);
      }


      std::vector<count_type> counts(mr.infos().nb_counts + 1);

      std::size_t n = mr.infos().nb_counts;

      mr.read<MAX_K, MAX_C>(kmer2, counts, n);

      bool f = true;

      while (kr.read<MAX_K, MAX_C>(kmer, count))
      {
        if (kmer < kmer2)
        {
          if (with_kmer)
          {
            kw->write<MAX_K, MAX_C>(kmer, count);
          }
          continue;
        }
        else if (kmer > kmer2)
        {
          if (f)
          {
            if (with_vector)
            {
              *vout << "0\n";
            }
            f = false;
          }
          while (mr.read<MAX_K, MAX_C>(kmer2, counts, n) && kmer > kmer2)
          {
            if (with_vector)
            {
              *vout << "0\n";
            }
          }

          if (kmer < kmer2)
          {
            if (with_kmer)
            {
              kw->write<MAX_K, MAX_C>(kmer, count);
            }
            f = true;
            continue;
          }
          else if (kmer == kmer2)
          {
            counts.back() = count;
            if (with_matrix)
            {
              mw->write<MAX_K, MAX_C>(kmer2, counts);
            }
            if (with_vector)
            {
              *vout << std::to_string(count) << '\n';
            }
          }
        }
        else
        {
          counts.back() = count;
          if (with_matrix)
          {
            mw->write<MAX_K, MAX_C>(kmer2, counts);
          }
          if (with_vector)
          {
            *vout << std::to_string(count) << '\n';
          }
        }

        f = false;
      }

      if (f && with_vector)
        *vout << "0\n";

      if (with_vector)
      {
        while (mr.read<MAX_K, MAX_C>(kmer2, counts, n))
        {
          *vout << "0\n";
        }
      }
    }

    void f_pa_matrix(const std::tuple<bool, bool, bool>& out_types)
    {
      const bool with_vector = std::get<0>(out_types);
      const bool with_matrix = std::get<1>(out_types);
      const bool with_kmer = std::get<2>(out_types);

      std::unique_ptr<std::ofstream> vout {nullptr};
      std::unique_ptr<PAMatrixWriter<8192>> mw {nullptr};
      std::unique_ptr<KmerWriter<8192>> kw {nullptr};

      KmerReader<8192> kr(m_kmers);
      Kmer<MAX_K> kmer; kmer.set_k(kr.infos().kmer_size);
      count_type count;

      PAMatrixReader<8192> mr(m_matrix);
      Kmer<MAX_K> kmer2; kmer2.set_k(mr.infos().kmer_size);

      if (with_vector)
      {
        vout = std::make_unique<std::ofstream>(m_vec, std::ios::out);
      }

      if (with_matrix)
      {
        mw = std::make_unique<PAMatrixWriter<8192>>(m_output,
                                                    mr.infos().kmer_size,
                                                    mr.infos().bits,
                                                    mr.infos().id,
                                                    mr.infos().partition,
                                                    m_cpr);
      }

      if (with_kmer)
      {
        kw = std::make_unique<KmerWriter<8192>>(m_koutput,
                                                kr.infos().kmer_size,
                                                kr.infos().count_slots,
                                                kr.infos().id,
                                                kr.infos().partition,
                                                m_cpr);
      }

      std::vector<uint8_t> bits(NBYTES(mr.infos().bits));

      mr.read<MAX_K>(kmer2, bits);

      bool f = true;

      while (kr.read<MAX_K, MAX_C>(kmer, count))
      {
        if (kmer < kmer2)
        {
          if (with_kmer)
          {
            kw->write<MAX_K, MAX_C>(kmer, count);
          }
          continue;
        }
        else if (kmer > kmer2)
        {
          if (f)
          {
            if (with_vector)
            {
              *vout << "0\n";
            }
            f = false;
          }
          while (mr.read<MAX_K>(kmer2, bits) && kmer > kmer2)
          {
            if (with_vector)
            {
              *vout << "0\n";
            }
          }

          if (kmer < kmer2)
          {
            if (with_kmer)
            {
              kw->write<MAX_K, MAX_C>(kmer, count);
            }
            f = true;
            continue;
          }
          else if (kmer == kmer2)
          {
            if (with_matrix)
            {
              mw->write<MAX_K>(kmer2, bits);
            }
            if (with_vector)
            {
              *vout << "1\n";
            }
          }
        }
        else if (kmer == kmer2)
        {
          if (with_matrix)
          {
            mw->write<MAX_K>(kmer2, bits);
          }
          if (with_vector)
          {
            *vout << "1\n";
          }
        }
        f = false;
      }

      if (f && with_vector)
        *vout << "0\n";

      if (with_vector)
      {
        while (mr.read<MAX_K>(kmer2, bits))
        {
          *vout << "0\n";
        }
      }
    }

  private:
    const std::string& m_matrix;
    const std::string& m_kmers;
    const std::string& m_output;
    const std::string& m_koutput;
    const std::string& m_vec;
    bool m_count;
    bool m_cpr;
    const std::tuple<bool, bool, bool>& m_out_types;
};

template<size_t MAX_K, size_t MAX_C>
class MatrixFilter
{
  using paths_t = std::vector<std::string>;

  public:
    MatrixFilter(const paths_t& matrices,
                 const paths_t& kmers,
                 const paths_t& outputs,
                 const paths_t& koutputs,
                 const paths_t& vecs,
                 bool cpr,
                 bool count,
                 std::size_t threads,
                 const std::tuple<bool, bool, bool>& out_types)
      : m_mpaths(matrices),
        m_kpaths(kmers),
        m_opaths(outputs),
        m_kopaths(koutputs),
        m_vopaths(vecs),
        m_cpr(cpr),
        m_count(count),
        m_threads(threads),
        m_out_types(out_types)
    {}

    void exec()
    {
      TaskPool pool(m_threads);
      for (std::size_t i = 0; i < m_mpaths.size(); i++)
      {
        task_t task = std::make_shared<FilterTask<MAX_K, MAX_C>>(m_mpaths[i], m_kpaths[i], m_opaths[i], m_kopaths[i], m_vopaths[i], m_cpr, m_count, m_out_types);

        pool.add_task(task);
      }

      pool.join_all();
    }

  private:
    const paths_t& m_mpaths;
    const paths_t& m_kpaths;
    const paths_t& m_opaths;
    const paths_t& m_kopaths;
    const paths_t& m_vopaths;

    bool m_cpr;
    bool m_count;
    std::size_t m_threads;
    const std::tuple<bool, bool, bool>& m_out_types;
};

} // namespace km
