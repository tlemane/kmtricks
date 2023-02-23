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
#include <algorithm>
#include <random>

#include <kmtricks/task.hpp>
#include <kmtricks/task_pool.hpp>
#include <kmtricks/cmd/all.hpp>
#include <kmtricks/gatb/gatb_utils.hpp>
#include <kmtricks/hash.hpp>
#include <kmtricks/progress.hpp>
#include <kmtricks/timer.hpp>
#include <indicators/progress_bar.hpp>
#include <indicators/cursor_control.hpp>


using namespace indicators;

namespace km {

template<size_t MAX_K, size_t MAX_C>
class TaskScheduler
{
public:
  TaskScheduler(all_options_t opt) : m_opt(opt), m_nb_samples(KmDir::get().m_fof.size())
  {
    if (spdlog::get_level() == spdlog::level::info)
      m_is_info = true;
    m_dyn.set_option(option::HideBarWhenComplete{false});
    init_progress();
  }

  ~TaskScheduler()
  {
    for (auto& p : m_progress) delete p;

    if (m_opt->hist)
    {
      for (auto& h : m_hists)
        HistWriter(KmDir::get().get_hist_path(KmDir::get().m_fof.get_id(h->idx())), *h, false);
    }
  }

  void init_progress()
  {
    m_progress.push_back(
      get_progress_bar("Configuration    ", 1, 50, Color::white, false));
    m_progress.push_back(
      get_progress_bar("Repartition      ", 1, 50, Color::white, false));
    m_progress.push_back(
      get_progress_bar("Compute SuperK   ", m_nb_samples, 50, Color::white, false));
  }

  void init_progress2(uint32_t nb_parts)
  {
    m_progress.push_back(
      get_progress_bar("Count partitions ", m_nb_samples * m_opt->restrict_to_list.size(), 50,
                        Color::white, false));
    m_progress.push_back(
      get_progress_bar("Merge partitions ", m_opt->restrict_to_list.size(), 50, Color::white, false));
    m_progress.push_back(
      get_progress_bar(
        "Format bloom     ", m_nb_samples, 50, Color::white, false));
  }

  void exec_config()
  {
    spdlog::info("Compute configuration...");
    IProperties* props = get_config_properties(m_opt->kmer_size,
                                               m_opt->minim_size,
                                               m_opt->minim_type,
                                               m_opt->repart_type,
                                               1,
                                               m_opt->nb_parts,
                                               m_opt->max_memory);
    ConfigTask<MAX_K> config_task(m_opt->fof, props, m_opt->bloom_size, m_opt->nb_parts);
    config_task.exec();
    Storage* config_storage = StorageFactory(STORAGE_FILE).load(KmDir::get().m_config_storage);
    LOCAL(config_storage);
    m_config.load(config_storage->getGroup("gatb"));
    KmDir::get().init_part(m_config._nb_partitions);

    m_hists.resize(m_nb_samples);
    for (size_t i=0; i<m_nb_samples; i++)
      m_hists[i] = m_opt->hist ? std::make_shared<KHist>(i, m_config._kmerSize, 1, 255) : nullptr;
  }

  void exec_repart()
  {
    spdlog::info("Compute minimizer repartition...");
    RepartTask<MAX_K> repart_task(m_opt->fof, m_opt->from);
    repart_task.exec(); repart_task.postprocess();
    m_opt->m_ab_min_vec.resize(KmDir::get().m_fof.size());
    m_hw = HashWindow(KmDir::get().m_hash_win);

    if (m_opt->restrict_to_list.empty())
    {
      if (m_opt->restrict_to != 1.0)
      {
        std::vector<uint32_t> ps;

        for (uint32_t i=0; i<m_config._nb_partitions; i++)
          ps.push_back(i);

        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(ps.begin(), ps.end(), g);
        size_t n = m_config._nb_partitions * m_opt->restrict_to;

        if (n == 0) n = 1;

        for (uint32_t i=0; i<n; i++)
        {
          m_opt->restrict_to_list.push_back(ps[i]);
        }
      }
      else
      {
        for (size_t p=0; p<m_config._nb_partitions; p++)
        {
          m_opt->restrict_to_list.push_back(p);
        }
      }
    }
    else
    {
      for (auto& p : m_opt->restrict_to_list)
      {
        if (p >= m_config._nb_partitions)
        {
          throw PipelineError(fmt::format("Ask to process part {} but nb_partitions is {}",
                                           p, m_config._nb_partitions));
        }
      }
    }
    init_progress2(m_config._nb_partitions);
  }

  void exec_superk()
  {
    if (m_is_info)
    {
      m_dyn.push_back(*m_progress[2]); m_dyn[0].set_progress(0);
    }

    TaskPool pool(m_opt->nb_threads);

    for (auto id : KmDir::get().m_fof)
    {
      task_t task = std::make_shared<SuperKTask<MAX_K>>(std::get<0>(id),
                                                        m_opt->lz4,
                                                        m_opt->restrict_to_list);
      if (m_is_info) task->set_callback([this](){ this->m_dyn[0].tick(); });

      spdlog::debug("[push] - SuperKTask - S={}", std::get<0>(id));
      pool.add_task(task);
    }
    pool.join_all();
    if (m_is_info) m_dyn[0].mark_as_completed();
  }

  void exec_count()
  {
    if (m_is_info) { m_dyn.push_back(*m_progress[3]); m_dyn[1].set_progress(0); }

    TaskPool pool(m_opt->nb_threads);

    for (auto id : KmDir::get().m_fof)
    {
      uint32_t a_min = std::get<2>(id) == 0 ? m_opt->c_ab_min : std::get<2>(id);
      uint32_t iid = KmDir::get().m_fof.get_i(std::get<0>(id));
      std::string sid = std::get<0>(id);
      sk_storage_t sk_storage = std::make_shared<SuperKStorageReader>(KmDir::get().get_superk_path(sid));
      parti_info_t pinfos = std::make_shared<PartiInfo<5>>(KmDir::get().get_superk_path(sid));
      for (auto& p : m_opt->restrict_to_list)
      {
        std::string path;
        task_t task = nullptr;
        if (m_opt->count_format == COUNT_FORMAT::KMER)
        {
          if (!m_opt->kff)
          {
            spdlog::debug("[push] - CountTask - S={}, P={}", sid, p);
            path = KmDir::get().get_count_part_path(
              sid, p, m_opt->lz4, KM_FILE::KMER);
            task = std::make_shared<CountTask<MAX_K, MAX_C, SuperKStorageReader>>(
              path, m_config, sk_storage, pinfos, p, iid, m_config._kmerSize,
              a_min, m_opt->lz4, get_hist_clone(m_hists[iid]), !m_opt->keep_tmp);
          }
          else if (m_opt->kff)
          {
            spdlog::debug("[push] - KffCountTask - S={}, P={}", sid, p);
            path = KmDir::get().get_count_part_path(
              sid, p, m_opt->lz4, KM_FILE::KFF);
            task = std::make_shared<KffCountTask<MAX_K, MAX_C, SuperKStorageReader>>(
              path, m_config, sk_storage, pinfos, p, iid,
              m_config._kmerSize, a_min, get_hist_clone(m_hists[iid]), !m_opt->keep_tmp);
          }
        }
        else
        {
          if (!m_opt->skip_merge)
          {
            spdlog::debug("[push] - HashCountTask - S={}, P={}", sid, p);
            path = KmDir::get().get_count_part_path(
              sid, p, m_opt->lz4, KM_FILE::HASH);
            task = std::make_shared<HashCountTask<MAX_K, MAX_C, SuperKStorageReader>>(
              path, m_config, sk_storage, pinfos, p, iid,
              m_hw.get_window_size_bits(), m_config._kmerSize, a_min, m_opt->lz4,
              get_hist_clone(m_hists[iid]), !m_opt->keep_tmp);
          }
          else
          {
            spdlog::debug("[push] - HashVecCountTask - S={}, P={}", sid, p);
            path = KmDir::get().get_count_part_path(
              sid, p, m_opt->lz4, KM_FILE::VECTOR);
            task = std::make_shared<HashVecCountTask<MAX_K, MAX_C, SuperKStorageReader>>(
              path, m_config, sk_storage, pinfos, p, iid,
              m_hw.get_window_size_bits(), m_config._kmerSize, a_min, m_opt->lz4,
              get_hist_clone(m_hists[iid]), !m_opt->keep_tmp);
          }
        }
        if (m_is_info) task->set_callback([this](){ this->m_dyn[1].tick(); });

        pool.add_task(task);
      }
    }
    pool.join_all();

    if (m_opt->hist)
    {
      for (auto& h : m_hists)
        h->merge_clones();
    }

    if (m_is_info) m_dyn[1].mark_as_completed();
  }

  void exec_superk_count()
  {
    if (m_is_info)
    {
      m_dyn.push_back(*m_progress[2]); m_dyn[0].set_progress(0);
      m_dyn.push_back(*m_progress[3]); m_dyn[1].set_progress(0);
    }
    TaskPool pool(m_opt->nb_threads);

    int max_running = std::floor(m_opt->nb_threads * m_opt->focus) > 0 ? m_opt->nb_threads * m_opt->focus : 1;

    for (auto id : KmDir::get().m_fof)
    {
      task_t task = std::make_shared<SuperKTask<MAX_K>>(std::get<0>(id),
                                                        m_opt->lz4,
                                                        m_opt->restrict_to_list);
      task->set_callback([this, id, &pool](){
        if (this->m_is_info)
          this->m_dyn[0].tick();
        uint32_t a_min = std::get<2>(id) == 0 ? this->m_opt->c_ab_min : std::get<2>(id);
        uint32_t iid = KmDir::get().m_fof.get_i(std::get<0>(id));
        std::string sid = std::get<0>(id);
        sk_storage_t sk_storage = std::make_shared<SuperKStorageReader>(KmDir::get().get_superk_path(sid));
        parti_info_t pinfos = std::make_shared<PartiInfo<5>>(KmDir::get().get_superk_path(sid));
        for (auto& p : this->m_opt->restrict_to_list)
        {
          std::string path;
          task_t task = nullptr;
          if (m_opt->count_format == COUNT_FORMAT::KMER)
          {
            if (!m_opt->kff)
            {
              spdlog::debug("[push] - CountTask - S={}, P={}", sid, p);
              path = KmDir::get().get_count_part_path(
                sid, p, this->m_opt->lz4, KM_FILE::KMER);
              task = std::make_shared<CountTask<MAX_K, MAX_C, SuperKStorageReader>>(
                path, this->m_config, sk_storage, pinfos, p, iid,
                this->m_config._kmerSize, a_min, m_opt->lz4, get_hist_clone(this->m_hists[iid]),
                !this->m_opt->keep_tmp);
            }
            else if (m_opt->kff)
            {
              spdlog::debug("[push] - KffCountTask - S={}, P={}", sid, p);
              path = KmDir::get().get_count_part_path(
                sid, p, this->m_opt->lz4, KM_FILE::KFF);
              task = std::make_shared<KffCountTask<MAX_K, MAX_C, SuperKStorageReader>>(
                path, this->m_config, sk_storage, pinfos, p, iid,
                this->m_config._kmerSize, a_min, get_hist_clone(this->m_hists[iid]), !this->m_opt->keep_tmp);
            }
          }
          else
          {
            if (!m_opt->skip_merge)
            {
              path = KmDir::get().get_count_part_path(
                sid, p, this->m_opt->lz4, KM_FILE::HASH);
              spdlog::debug("[push] - HashCountTask - S={}, P={}", sid, p);
              task = std::make_shared<HashCountTask<MAX_K, MAX_C, SuperKStorageReader>>(
                path, m_config, sk_storage, pinfos, p, iid,
                m_hw.get_window_size_bits(), m_config._kmerSize, a_min, m_opt->lz4,
                get_hist_clone(this->m_hists[iid]), !this->m_opt->keep_tmp);
            }
            else
            {
              spdlog::debug("[push] - HashVecCountTask - S={}, P={}", sid, p);
              path = KmDir::get().get_count_part_path(
                sid, p, false, KM_FILE::VECTOR);
              task = std::make_shared<HashVecCountTask<MAX_K, MAX_C, SuperKStorageReader>>(
                path, this->m_config, sk_storage, pinfos, p, iid,
                this->m_hw.get_window_size_bits(), this->m_config._kmerSize, a_min, false,
                get_hist_clone(this->m_hists[iid]), !this->m_opt->keep_tmp);
            }
          }
          if (m_is_info)
          {
            ProgressBar* ptr = &this->m_dyn[1];
            task->set_callback([ptr](){ ptr->tick(); });
          }
          pool.add_task(task);
        }
      });

      while (superk_in() >= max_running)
      {
        if (max_running == m_opt->nb_threads)
        {
          max_running /= 2;
          if (!(max_running > 0)) max_running = 1;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
      }
      task->set_level(5);
      m_superk.push_back(task);
      spdlog::debug("[push] - SuperKTask - S={}", std::get<0>(id));
      pool.add_task(task);
    }
    while (superk_finish() != m_nb_samples)
    {
      std::this_thread::sleep_for(std::chrono::milliseconds(20));
    }
    pool.join_all();

    if (m_opt->hist)
    {
      for (auto& h : m_hists)
        h->merge_clones();
    }

    if (m_is_info)
      m_dyn[0].mark_as_completed();
  }

  size_t superk_finish()
  {
    size_t count = 0;
    for (auto& e : m_superk)
    {
      if (e->finish()) count++;
    }
    return count;
  }

  size_t superk_running()
  {
    size_t count = 0;
    for (auto& e : m_superk)
    {
      if (e->running()) count++;
    }
    return count;
  }


  int superk_in()
  {
    int count = 0;
    for (auto& e : m_superk)
    {
      if (e->in_queue()) count++;
    }
    return count;
  }

  void exec_merge()
  {
    if (m_is_info)
    {
      m_dyn.push_back(*m_progress[4]); m_dyn[2].set_progress(0);
    }

    if (m_opt->m_ab_float)
    {
      m_opt->m_ab_min_vec = compute_merge_thresholds(m_hists, m_opt->m_ab_min_f,
                                                     KmDir::get().get_merge_th_path());
    }
    TaskPool pool(m_opt->nb_threads);
    for (auto& p : m_opt->restrict_to_list)
    {
      task_t task = nullptr;
      if (m_opt->count_format == COUNT_FORMAT::KMER)
      {
        spdlog::debug("[push] - KmerMergeTask - P={}", p);
        task = std::make_shared<KmerMergeTask<MAX_K, MAX_C>>(
          p, m_opt->m_ab_min_vec, m_config._kmerSize, m_opt->r_min, m_opt->save_if,
          m_opt->lz4, m_opt->mode, m_opt->format, !m_opt->keep_tmp);
      }
      else if (m_opt->count_format == COUNT_FORMAT::HASH)
      {
        spdlog::debug("[push] - HashMergeTask - P={}", p);
        task = std::make_shared<HashMergeTask<MAX_C>>(
          p, m_opt->m_ab_min_vec, m_opt->r_min, m_opt->save_if, m_opt->lz4, m_opt->mode,
          m_opt->format, m_hw, !m_opt->keep_tmp, m_opt->bwidth);
      }
      if (m_is_info) task->set_callback([this](){ this->m_dyn[2].tick(); });
      pool.add_task(task);
    }
    pool.join_all();
    if (m_is_info)
      m_dyn[2].mark_as_completed();
  }

  void exec_format()
  {
    TaskPool pool(m_opt->nb_threads);
    if (m_is_info)
    {
      m_dyn.push_back(*m_progress[5]);
      if (m_opt->skip_merge) m_dyn[2].set_progress(0);
      else m_dyn[3].set_progress(0);
    }

    if (m_opt->skip_merge)
    {
      for (auto id : KmDir::get().m_fof)
      {
        spdlog::debug("[push] - FormatVectorTask - S={}", std::get<0>(id));
        task_t task = std::make_shared<FormatVectorTask>(
          std::get<0>(id), m_opt->out_format, m_hw.bloom_size(), m_config._nb_partitions, false, m_config._kmerSize, !m_opt->keep_tmp);
        if (m_is_info)
          task->set_callback([this](){ this->m_dyn[2].tick(); });
        pool.add_task(task);
      }
      pool.join_all();

      if (m_is_info)
        m_dyn[2].mark_as_completed();
    }
    else
    {
      std::vector<vmr_t<8192>> matrix_files; matrix_files.reserve(m_config._nb_partitions);
      std::vector<int> fds; fds.reserve(m_config._nb_partitions);
      std::vector<std::mutex> mutex(m_config._nb_partitions);
      for (size_t p=0; p<m_config._nb_partitions; p++)
      {
        fds.push_back(open(KmDir::get().get_matrix_path(p, MODE::BFT, FORMAT::BIN, COUNT_FORMAT::HASH, false).c_str(), O_RDONLY));
      }

      for (auto id : KmDir::get().m_fof)
      {
        spdlog::debug("[push] - FormatTask - S={}", std::get<0>(id));
        std::string sid = std::get<0>(id);
        uint32_t file_id = KmDir::get().m_fof.get_i(sid);
        task_t task = std::make_shared<FormatTask>(
          fds, mutex, m_opt->out_format, m_hw.bloom_size(), file_id, m_config._nb_partitions,
          m_config._kmerSize, !m_opt->keep_tmp);
        if (m_is_info)
          task->set_callback([this](){ this->m_dyn[3].tick(); });
        pool.add_task(task);
      }
      pool.join_all();

      if (!m_opt->keep_tmp)
      {
        for (auto s: KmDir::get().get_matrix_paths(m_opt->restrict_to_list.size(), MODE::BFT, FORMAT::BIN,
                                                    COUNT_FORMAT::HASH, false))
        {
          Eraser::get().erase(s);
        }
      }

      if (m_is_info)
        m_dyn[3].mark_as_completed();
    }
  }

  void execute()
  {
    Timer whole_time;

    exec_config();
    exec_repart();

    if (m_opt->until == COMMAND::REPART)
      goto end;

    if (m_opt->until == COMMAND::SUPERK)
    {
      exec_superk();
      goto end;
    }

    exec_superk_count();

    if (m_opt->until == COMMAND::COUNT)
      goto end;

    if (!m_opt->skip_merge && !m_opt->kff)
    {
      exec_merge();

      if (m_opt->until == COMMAND::MERGE)
        goto end;
    }

    if (m_opt->mode == MODE::BFT)
      exec_format();

    end:
      spdlog::info("Done in {} - Peak RSS -> {:.2f} MB.",
                   whole_time.formatted(),
                   get_peak_rss() * 0.0009765625);

    std::ofstream out_infos(KmDir::get().m_run_infos, std::ios::out);
    check_fstream_good(KmDir::get().m_run_infos, out_infos);
    out_infos << "Time: " << std::to_string(whole_time.elapsed<std::chrono::seconds>().count());
    out_infos << " seconds" << "\n";
    out_infos << "Memory: " << std::to_string(get_peak_rss() * 0.0009765625) << "MB" << std::endl;
    Eraser::get().join();
    return;
  }

public:
  all_options_t m_opt;
  Configuration m_config;
  std::vector<task_t> m_superk;
  std::vector<task_t> m_counts;
  std::vector<hist_t> m_hists;
  TaskPool* m_pool {nullptr};
  std::condition_variable m_cv;
  size_t m_nb_samples;
  HashWindow m_hw;
  bool m_is_info {false};

  std::vector<ProgressBar*> m_progress;
  DynamicProgress<ProgressBar> m_dyn;
};

};
