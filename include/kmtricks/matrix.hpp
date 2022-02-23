#pragma once

#include <vector>
#include <string>
#include <filesystem>

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
               bool cpr, bool count)
      : ITask(0),
        m_matrix(matrix),
        m_kmers(kmers),
        m_output(output),
        m_koutput(koutput),
        m_vec(vec),
        m_cpr(cpr),
        m_count(count)
    {
    }

    void exec()
    {
      if (m_count)
        f_count_matrix();
      else
        f_pa_matrix();
    }

    void preprocess() {}
    void postprocess()
    {
      fs::remove(m_kmers);
    }

  private:

    void f_count_matrix()
    {
      std::ofstream vout(m_vec, std::ios::out);

      KmerReader<8192> kr(m_kmers);
      Kmer<MAX_K> kmer; kmer.set_k(kr.infos().kmer_size);
      count_type count;

      MatrixReader<8192> mr(m_matrix);
      Kmer<MAX_K> kmer2; kmer2.set_k(mr.infos().kmer_size);

      MatrixWriter<8192> mw(m_output,
                              mr.infos().kmer_size,
                              mr.infos().count_slots,
                              mr.infos().nb_counts + 1,
                              mr.infos().id,
                              mr.infos().partition,
                              m_cpr);

      KmerWriter<8192> kw(m_koutput,
                          kr.infos().kmer_size,
                          kr.infos().count_slots,
                          kr.infos().id,
                          kr.infos().partition,
                          m_cpr);

      std::vector<count_type> counts(mr.infos().nb_counts + 1);

      std::size_t n = mr.infos().nb_counts;

      mr.read<MAX_K, MAX_C>(kmer2, counts, n);

      bool f = true;

      while (kr.read<MAX_K, MAX_C>(kmer, count))
      {
        if (kmer < kmer2)
        {
          kw.write<MAX_K, MAX_C>(kmer, count);
        }
        else if (kmer > kmer2)
        {
          if (f)
          {
            vout << "0\n";
            f = false;
          }
          while (mr.read<MAX_K, MAX_C>(kmer2, counts, n) && kmer > kmer2)
          {
            vout << "0\n";
          }

          if (kmer < kmer2)
          {
            kw.write<MAX_K, MAX_C>(kmer, count);
            vout << "0\n";
          }
          else
          {
            counts.back() = count;
            mw.write<MAX_K, MAX_C>(kmer2, counts);
            vout << std::to_string(count) << '\n';
          }
        }
        else
        {
          counts.back() = count;
          mw.write<MAX_K, MAX_C>(kmer2, counts);
          vout << std::to_string(count) << '\n';
        }

        f = false;
      }

      while (mr.read<MAX_K, MAX_C>(kmer2, counts, n))
      {
        vout << "0\n";
      }
    }

    void f_pa_matrix()
    {
      std::ofstream vout(m_vec, std::ios::out);

      KmerReader<8192> kr(m_kmers);
      Kmer<MAX_K> kmer; kmer.set_k(kr.infos().kmer_size);
      count_type count;

      PAMatrixReader<8192> mr(m_matrix);
      Kmer<MAX_K> kmer2; kmer2.set_k(mr.infos().kmer_size);

      PAMatrixWriter<8192> mw(m_output,
                              mr.infos().kmer_size,
                              mr.infos().bits,
                              mr.infos().id,
                              mr.infos().partition,
                              m_cpr);

      KmerWriter<8192> kw(m_koutput,
                          kr.infos().kmer_size,
                          kr.infos().count_slots,
                          kr.infos().id,
                          kr.infos().partition,
                          m_cpr);

      std::vector<uint8_t> bits(NBYTES(mr.infos().bits));

      mr.read<MAX_K>(kmer2, bits);

      bool f = true;

      while (kr.read<MAX_K, MAX_C>(kmer, count))
      {
        if (kmer < kmer2)
        {
          kw.write<MAX_K, MAX_C>(kmer, count);
          continue;
        }
        else if (kmer > kmer2)
        {
          if (f)
          {
            vout << "0\n";
            //vout << kmer2.to_string() << " 0\n";
            f = false;
          }
          while (mr.read<MAX_K>(kmer2, bits) && kmer > kmer2)
          {
            vout << "0\n";
            //vout << kmer2.to_string() << " 0\n";
          }

          if (kmer < kmer2)
          {
            kw.write<MAX_K, MAX_C>(kmer, count);
            f = true;
            continue;
          }
          else if (kmer == kmer2)
          {
            mw.write<MAX_K>(kmer2, bits);
            vout << "1\n";
            //vout << kmer2.to_string() << " 1\n";
          }
        }
        else if (kmer == kmer2)
        {
          mw.write<MAX_K>(kmer2, bits);
          vout << "1\n";
          //vout << kmer2.to_string() << " 1\n";
        }
        f = false;
      }

      if (f)
        vout << "0\n";
        //vout << kmer2.to_string() << " 0\n";
      while (mr.read<MAX_K>(kmer2, bits))
      {
        vout << "0\n";
        //vout << kmer2.to_string() << " 0\n";
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
                 std::size_t threads)
      : m_mpaths(matrices),
        m_kpaths(kmers),
        m_opaths(outputs),
        m_kopaths(koutputs),
        m_vopaths(vecs),
        m_cpr(cpr),
        m_count(count),
        m_threads(threads)
    {}

    void exec()
    {
      TaskPool pool(m_threads);
      for (std::size_t i = 0; i < m_mpaths.size(); i++)
      {
        task_t task = std::make_shared<FilterTask<MAX_K, MAX_C>>(m_mpaths[i], m_kpaths[i], m_opaths[i], m_kopaths[i], m_vopaths[i], m_cpr, m_count);

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
};

} // namespace km
