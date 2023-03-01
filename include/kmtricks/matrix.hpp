#pragma once

#include <vector>
#include <string>
#include <filesystem>
#include <tuple>
#include <memory>
#include <sstream>

#include <kmtricks/utils.hpp>
#include <kmtricks/io/pa_matrix_file.hpp>
#include <kmtricks/io/kmer_file.hpp>
#include <kmtricks/itask.hpp>
#include <kmtricks/task_pool.hpp>
#include <kmtricks/hash.hpp>
#include <kmtricks/repartition.hpp>
#include <kmtricks/io/fof.hpp>

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
          //
          else
          {
            if (with_kmer)
            {
              kw->write<MAX_K, MAX_C>(kmer, count);
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
          else
          {
            if (with_kmer)
            {
              kw->write<MAX_K, MAX_C>(kmer, count);
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

template<std::size_t MAX_K, std::size_t MAX_C>
class MatrixMergeTask;

template<std::size_t MAX_K, std::size_t MAX_C>
class MatrixMerger
{
  enum class mmode : std::uint8_t {
    kmer, hash
  };

  static constexpr mmode mode = (MAX_K == 1) ? mmode::hash : mmode::kmer;

  using kmer_type = std::conditional_t<
    mode == mmode::kmer, Kmer<MAX_K>, std::uint64_t
  >;

  using count_type = std::conditional_t<
    MAX_C == 1, std::uint8_t, typename selectC<MAX_C>::type
  >;

  using data_type = std::vector<count_type>;

  static constexpr std::size_t buf_size = 8192;

  using input_stream_type = std::unique_ptr<
    std::conditional_t<
      mode == mmode::kmer,
      std::conditional_t<
        MAX_C == 1, PAMatrixReader<buf_size>, MatrixReader<buf_size>>,
      std::conditional_t<
        MAX_C == 1, PAHashMatrixReader<buf_size>, MatrixHashReader<buf_size>>
    >
  >;

  using output_stream_type = std::unique_ptr<
    std::conditional_t<
      mode == mmode::kmer,
      std::conditional_t<
        MAX_C == 1, PAMatrixWriter<buf_size>, MatrixWriter<buf_size>>,
      std::conditional_t<
        MAX_C == 1, PAHashMatrixWriter<buf_size>, MatrixHashWriter<buf_size>>
    >
  >;

  public:
    class PartitionMerger
    {
      struct element
      {
        kmer_type value;
        data_type data;
        std::size_t pos;
        std::size_t n;
        input_stream_type stream;
        bool is_set {false};


        element(const std::string& path, std::size_t p)
          : pos(p)
        {
          bool is_kmer = (
            fs::path(path).filename().string().find("kmer") != std::string::npos
          );
          if constexpr(mode == mmode::kmer && MAX_C != 1)
            stream = std::make_unique<typename input_stream_type::element_type>(path, is_kmer);
          else
            stream = std::make_unique<typename input_stream_type::element_type>(path);

          if constexpr(MAX_C != 1)
          {
            n = stream->infos().nb_counts;
            data.resize(n);
          }
          else
          {
            n = stream->infos().bits;
            data.resize(stream->infos().bytes);
          }

          if constexpr(mode == mmode::kmer)
            value.set_k(stream->infos().kmer_size);

          load();
        }

        void load()
        {
          if constexpr(mode == mmode::kmer)
          {
            if constexpr(MAX_C != 1)
              is_set = stream->template read<MAX_K, MAX_C>(value, data);
            else
              is_set = stream->template read<MAX_K>(value, data);
          }
          else
          {
            if constexpr(MAX_C != 1)
              is_set = stream->template read<MAX_C>(value, data);
            else
              is_set = stream->read(value, data);
          }
        }

        explicit operator bool() const { return is_set; }
      };

      using element_type = std::unique_ptr<element>;

      struct element_cmp
      {
        bool operator()(element* lhs, element* rhs) const {
          return lhs->value > rhs->value;
        }
      };

      using queue_type = std::priority_queue<
        element*, std::vector<element*>, element_cmp
      >;

      public:

        PartitionMerger() : m_init(false) {}

        PartitionMerger(const std::vector<std::string>& paths)
        {
          init(paths);
        }

        void init(const std::vector<std::string>& paths)
        {
          if (m_init)
            return;

          std::size_t pos = 0;
          for (auto& p : paths)
          {
            m_elements.push_back(std::make_unique<element>(p, pos));
            m_queue.push(m_elements.back().get());
            pos += m_elements.back()->n;
          }
          if constexpr (MAX_C != 1)
            m_current_data.resize(pos);
          else
            m_current_data.resize(NBYTES(pos));
          m_init = true;
        }

        bool next()
        {
          if (m_queue.empty())
            return false;

          std::fill(m_current_data.begin(), m_current_data.end(), 0);
          auto elem = m_queue.top();
          m_current_kmer = elem->value;

          if constexpr(MAX_C != 1)
            std::copy(elem->data.begin(), elem->data.end(), m_current_data.begin() + elem->pos);
          else
            copy_pa_vec(elem->pos, elem->n, elem->data);

          m_queue.pop();

          elem->load();
          if (*elem)
            m_queue.push(elem);

          if (m_queue.empty())
            return false;

          for (elem = m_queue.top(); elem->value == m_current_kmer; elem = m_queue.top())
          {
            if constexpr(MAX_C != 1)
              std::copy(elem->data.begin(), elem->data.end(), m_current_data.begin() + elem->pos);
            else
              copy_pa_vec(elem->pos, elem->n, elem->data);

            m_queue.pop();
            elem->load();

            if (*elem)
              m_queue.push(elem);

            if (m_queue.empty())
              return true;
          }

          return true;
        }

        const kmer_type& current_kmer() const { return m_current_kmer; }
        const data_type& current_data() const { return m_current_data; }

        void write(const std::string& path, bool cpr)
        {
          if constexpr(mode == mmode::kmer)
          {
            if constexpr(MAX_C == 1)
              write_k_p(path, cpr);
            else
              write_k_c(path, cpr);
          }
          else
          {
            if constexpr(MAX_C == 1)
              write_h_p(path, cpr);
            else
              write_h_c(path, cpr);
          }
        }

      private:

        // TODO This is only temporary for testing, we have to use bit packing
        void copy_pa_vec(std::size_t start, std::size_t n, const data_type& data)
        {
          std::size_t j = 0;
          for (std::size_t i = start; i < start + n; ++i)
          {
            if (BITCHECK(data, j))
              BITSET(m_current_data, i);
            ++j;
          }
        }

        std::size_t get_ns() const
        {
          std::size_t n = 0;
          for (auto& e : m_elements)
            n += e->n;
          return n;
        }

        void write_k_c(const std::string& path, bool cpr)
        {
          auto& i = m_elements.back()->stream->infos();
          auto out = std::make_unique<typename output_stream_type::element_type>(
            path, i.kmer_size, i.count_slots, get_ns(), i.id, i.partition, cpr
          );

          while (next())
            out->template write<MAX_K, MAX_C>(m_current_kmer, m_current_data);
        }

        void write_k_p(const std::string& path, bool cpr)
        {
          auto& i = m_elements.back()->stream->infos();
          auto out = std::make_unique<typename output_stream_type::element_type>(
            path, i.kmer_size, get_ns(), i.id, i.partition, cpr
          );

          while (next())
            out->template write<MAX_K>(m_current_kmer, m_current_data);
        }

        void write_h_c(const std::string& path, bool cpr)
        {
          auto& i = m_elements.back()->stream->infos();
          auto out = std::make_unique<typename output_stream_type::element_type>(
            path, i.count_slots, get_ns(), i.id, i.partition, cpr
          );

          while (next())
            out->template write<MAX_C>(m_current_kmer, m_current_data);
        }

        void write_h_p(const std::string& path, bool cpr)
        {
          auto& i = m_elements.back()->stream->infos();
          auto out = std::make_unique<typename output_stream_type::element_type>(
            path, get_ns(), i.id, i.partition, cpr
          );

          while (next())
            out->write(m_current_kmer, m_current_data);
        }


      private:
        kmer_type m_current_kmer;
        data_type m_current_data;
        queue_type m_queue;
        std::vector<element_type> m_elements;
        bool m_init {false};
    };
  public:

    MatrixMerger(const std::vector<std::string>& runs, const std::string& output, bool cpr)
      : m_runs(runs), m_output(output), m_cpr(cpr)
    {
      sanity_check();
      copy_km_dir();
      init_nb_part();
      cat_fof();
    }

    std::shared_ptr<MatrixMergeTask<MAX_K, MAX_C>> make_task(std::size_t p) const
    {
      std::vector<std::string> paths = paths_from_runs(p);
      return std::make_shared<MatrixMergeTask<MAX_K, MAX_C>>(
        PartitionMerger(paths), output_path(p), m_cpr
      );
    }

    void exec(TaskPool& pool)
    {
      for (std::size_t p = 0; p < nb_parts(); ++p)
        pool.add_task(make_task(p));
      pool.join_all();
    }

    std::vector<std::string> get_merge_paths(std::size_t p) const
    {
      return paths_from_runs(p);
    }

    std::size_t nb_parts() const
    {
      return m_nb_parts;
    }
  private:

    void sanity_check() const
    {
      Repartition r1(m_runs[0] + "/repartition_gatb/repartition.minimRepart");
      auto& table = r1.table();

      if (!fs::exists(m_runs[0] + "/repartition_gatb/repartition.minimRepart"))
        throw InputError(m_runs[0] + ": not a kmtricks directory.");

      for (std::size_t i = 1; i < m_runs.size(); ++i)
      {
        Repartition r2(m_runs[i] + "/repartition_gatb/repartition.minimRepart");
        auto& table2 = r2.table();

        if (!std::equal(table.begin(), table.end(), table2.begin(), table2.end()))
          throw InputError(m_runs[0] + " and " + m_runs[i] + " are not mergeable.") ;
      }
    }

    void copy_km_dir() const
    {
      auto& p = m_runs[0];
      fs::create_directory(m_output);
      fs::create_directory(m_output + "/matrices");
      fs::create_directory(m_output + "/repartition_gatb");
      fs::create_directory(m_output + "/config_gatb");

      fs::copy(p + "/hash.info", m_output, fs::copy_options::recursive);
      fs::copy(p + "/config_gatb", m_output + "/config_gatb", fs::copy_options::recursive);
      fs::copy(p + "/repartition_gatb", m_output + "/repartition_gatb", fs::copy_options::recursive);
      fs::copy(p + "/options.txt", m_output, fs::copy_options::recursive);
    }

    std::vector<std::string> paths_from_runs(std::size_t p) const
    {
      std::vector<std::string> paths;

      for (auto& r : m_runs)
      {
        auto is_kmer = is_kmer_run(r);
        if (is_kmer)
        {
          auto kp = get_k_paths(r, p);
          paths.insert(
            paths.end(), std::make_move_iterator(kp.begin()), std::make_move_iterator(kp.end()));

        }
        else
        {
          paths.push_back(get_mat_path(r, p));
        }
      }
      return paths;
    }

    std::vector<std::string> get_k_paths(const std::string& run, std::size_t p) const
    {
      std::stringstream mp;
      mp << run << "/counts/partition_" << std::to_string(p);

      const fs::path cpath(mp.str());
      std::vector<std::string> paths;
      for (auto& entry : fs::directory_iterator(cpath))
        paths.push_back(entry.path().string());
      return paths;
    }

    std::string get_mat_path(const std::string& run, std::size_t p) const
    {
      std::stringstream mp;
      mp << run << "/matrices/";

      const fs::path cpath(mp.str());
      std::vector<std::string> paths;
      for (auto& entry : fs::directory_iterator(cpath))
        paths.push_back(entry.path().string());

      std::sort(std::begin(paths), std::end(paths));
      return paths[p];
    }

    bool is_kmer_run(const std::string& run) const
    {
      std::stringstream ss;
      ss << run << "/counts/partition_0/";
      return !fs::is_empty(fs::path(ss.str()));
    }

    std::string output_path(std::size_t p) const
    {
      std::string path = m_output;
      path.append("/matrices/matrix_" + std::to_string(p));

      if (mode == mmode::kmer && MAX_C != 1)
        path.append(".count");
      else if (mode == mmode::kmer && MAX_C == 1)
        path.append(".pa");
      else if (mode == mmode::hash && MAX_C != 1)
        path.append(".count_hash");
      else if (mode == mmode::hash && MAX_C == 1)
        path.append(".pa_hash");

      if (m_cpr)
        path.append(".lz4");

      return path;
    }

    void cat_fof() const
    {
      {
        std::ofstream out(m_output + "/kmtricks.fof", std::ios::out);

        for (auto& r : m_runs)
        {
          std::ifstream inf(r + "/kmtricks.fof", std::ios::in);

          for (std::string line; std::getline(inf, line);)
          {
            if (!line.empty())
              out << line << '\n';
          }
        }
      }

      try {
        Fof(m_output + "/kmtricks.fof");
      } catch (...) {
        cat_fof_and_rename();
      }
    }

    void cat_fof_and_rename() const
    {
      {
        std::ofstream out(m_output + "/kmtricks.fof", std::ios::out);

        std::size_t nr = 0;
        for (auto& r : m_runs)
        {
          std::ifstream inf(r + "/kmtricks.fof", std::ios::in);

          for (std::string line; std::getline(inf, line);)
          {
            if (!line.empty())
            {
              std::stringstream ss;
              auto ls = bc::utils::split(line, ':');
              auto id = bc::utils::trim(ls[0]);
              out << id << "_" << std::to_string(nr) << ": " << ls[1] << '\n';
            }
          }
          ++nr;
        }
      }
    }

    void init_nb_part()
    {
      std::ifstream inf(m_runs[0] + "/hash.info", std::ios::in | std::ios::binary);
      inf.seekg(8);
      inf.read(reinterpret_cast<char*>(&m_nb_parts), sizeof(m_nb_parts));
    }

  private:
    std::vector<std::string> m_runs;
    std::string m_output;
    bool m_cpr;
    std::size_t m_nb_parts {0};
};

template<std::size_t MAX_K, std::size_t MAX_C>
class MatrixMergeTask : public ITask
{
  using partition_merge_type = typename MatrixMerger<MAX_K, MAX_C>::PartitionMerger;

  public:
    MatrixMergeTask(partition_merge_type&& pm, const std::string& output, bool cpr)
      : ITask(0), m_pm(std::move(pm)), m_output(output), m_cpr(cpr)
    {

    }

    void exec() override
    {
      m_pm.write(m_output, m_cpr);
    }

    void preprocess() override {}
    void postprocess() override {}


  private:
   partition_merge_type m_pm;
   std::string m_output;
   bool m_cpr;
};

} // namespace km
