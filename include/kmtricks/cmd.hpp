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

#include <kmtricks/cmd/cmd_common.hpp>
#include <kmtricks/cmd/all.hpp>
#include <kmtricks/cmd/repart.hpp>
#include <kmtricks/cmd/superk.hpp>
#include <kmtricks/cmd/count.hpp>
#include <kmtricks/cmd/merge.hpp>
#include <kmtricks/cmd/format.hpp>
#include <kmtricks/cmd/dump.hpp>
#include <kmtricks/cmd/infos.hpp>
#include <kmtricks/cmd/aggregate.hpp>
#include <kmtricks/cmd/filter.hpp>
#include <kmtricks/cmd/index.hpp>
#include <kmtricks/cmd/query.hpp>
#include <kmtricks/cmd/combine.hpp>

#include <kmtricks/io.hpp>
#include <kmtricks/utils.hpp>
#include <kmtricks/task.hpp>
#include <kmtricks/gatb/gatb_utils.hpp>
#include <kmtricks/kmdir.hpp>
#include <kmtricks/task_pool.hpp>
#include <kmtricks/task_scheduler.hpp>
#include <kmtricks/progress.hpp>
#include <kmtricks/signals.hpp>
#include <kmtricks/matrix.hpp>

#ifdef WITH_PLUGIN
#include <kmtricks/plugin_manager.hpp>
#include <kmtricks/plugin.hpp>
#endif

#ifdef WITH_HOWDE
#include <cmd_cluster.h>
#include <cmd_build_sbt.h>
#include <cmd_query.h>
#endif

namespace km {

template<size_t MAX_K>
struct main_all
{
  void operator()(km_options_t options)
  {
    spdlog::info("Run with {} implementation", Kmer<MAX_K>::name());
    all_options_t opt = std::static_pointer_cast<struct all_options>(options);
    spdlog::debug(opt->display());
    opt->sanity_check();
    KmDir::get().init(opt->dir, opt->fof, true);
    opt->dump(KmDir::get().m_options);

#ifdef WITH_PLUGIN
    if (opt->use_plugin)
      PluginManager<IMergePlugin>::get().init(opt->plugin, opt->plugin_config, MAX_K);
#endif

    TaskScheduler<MAX_K, DMAX_C> scheduler(opt);
    scheduler.execute();
  }
};

template<size_t MAX_K>
struct main_repart
{
  void operator()(km_options_t options)
  {
    spdlog::info("Run with {} implementation", Kmer<MAX_K>::name());
    repart_options_t opt = std::static_pointer_cast<struct repart_options>(options);
    spdlog::debug(opt->display());
    KmDir::get().init(opt->dir, opt->fof, true);
    IProperties* props = get_config_properties(opt->kmer_size,
                                          opt->minim_size,
                                          opt->minim_type,
                                          opt->repart_type,
                                          1,
                                          opt->nb_parts);
    ConfigTask<MAX_K> config_task(opt->fof, props, opt->bloom_size, opt->nb_parts);
    config_task.exec();
    RepartTask<MAX_K> repart_task(opt->fof); repart_task.exec(); repart_task.postprocess();

    Storage* config_storage = StorageFactory(STORAGE_FILE).load(KmDir::get().m_config_storage);
    LOCAL(config_storage);

    Configuration config = Configuration();
    config.load(config_storage->getGroup("gatb"));
    KmDir::get().init_part(config._nb_partitions);
  }
};

template<size_t MAX_K>
struct main_superk
{
  void operator()(km_options_t options)
  {
    spdlog::info("Run with {} implementation", Kmer<MAX_K>::name());
    superk_options_t opt = std::static_pointer_cast<struct superk_options>(options);
    spdlog::debug(opt->display());
    KmDir::get().init(opt->dir, "", false);

    Storage* config_storage = StorageFactory(STORAGE_FILE).load(KmDir::get().m_config_storage);
    LOCAL(config_storage);
    Configuration config = Configuration();
    config.load(config_storage->getGroup("gatb"));

    if (opt->restrict_to_list.empty())
    {
      for (uint32_t p=0; p<config._nb_partitions; p++)
      {
        opt->restrict_to_list.push_back(p);
      }
    }
    else
    {
      for (auto& p : opt->restrict_to_list)
        if (p >= config._nb_partitions)
          throw ConfigError(fmt::format("Ask to process partition {} but nb_partitions is {}",
                                        p, config._nb_partitions));
    }

    SuperKTask<MAX_K> superk_task(opt->id, opt->lz4, opt->restrict_to_list); superk_task.exec();
  }
};

template<size_t MAX_K>
struct main_count
{
  void operator()(km_options_t options)
  {
    spdlog::info("Run with {} implementation", Kmer<MAX_K>::name());
    count_options_t opt = std::static_pointer_cast<struct count_options>(options);
    spdlog::debug(opt->display());
    KmDir::get().init(opt->dir, "", false);

    Storage* config_storage = StorageFactory(STORAGE_FILE).load(KmDir::get().m_config_storage);
    LOCAL(config_storage);
    Configuration config = Configuration();
    config.load(config_storage->getGroup("gatb"));

    sk_storage_t superk_storage = std::make_shared<SuperKStorageReader>(
      KmDir::get().get_superk_path(opt->id));
    parti_info_t pinfo = std::make_shared<PartiInfo<5>>(KmDir::get().get_superk_path(opt->id));

    TaskPool pool(opt->nb_threads);
    HashWindow hw(KmDir::get().m_hash_win);

    hist_t hist = opt->hist ? std::make_shared<KHist>(KmDir::get().m_fof.get_i(opt->id),
                                          config._kmerSize, 1, 255) : nullptr;

    for (size_t i=0; i<config._nb_partitions; i++)
    {
      if (opt->partition_id != -1)
        if (static_cast<size_t>(opt->partition_id) != i)
          continue;
      if (opt->format == "kmer" || opt->format == "kff")
      {
        std::string path = KmDir::get().get_count_part_path(
          opt->id, i, opt->lz4, opt->format == "kmer" ? KM_FILE::KMER : KM_FILE::KFF);

        if (opt->format == "kmer")
        {
          spdlog::debug("[push] - CountTask - S={}, P={}", opt->id, i);
          pool.add_task(std::make_shared<CountTask<MAX_K, DMAX_C, SuperKStorageReader>>(
            path, config, superk_storage, pinfo, i, KmDir::get().m_fof.get_i(opt->id),
            config._kmerSize, opt->c_ab_min, opt->lz4, get_hist_clone(hist), opt->clear));
        }
        else if (opt->format == "kff")
        {
          spdlog::debug("[push] - KffCountTask - S={}, P={}", opt->id, i);
          pool.add_task(std::make_shared<KffCountTask<MAX_K, DMAX_C, SuperKStorageReader>>(
            path, config, superk_storage, pinfo, i, KmDir::get().m_fof.get_i(opt->id),
            config._kmerSize, opt->c_ab_min, get_hist_clone(hist), opt->clear));
        }
      }
      else if (opt->format == "hash" || opt->format == "vector")
      {
        std::string path = KmDir::get().get_count_part_path(
          opt->id, i, opt->lz4, opt->format == "hash" ? KM_FILE::HASH : KM_FILE::VECTOR);

        if (opt->format == "hash")
        {
          spdlog::debug("[push] - HashCountTask - S={}, P={}", opt->id, i);
          pool.add_task(std::make_shared<HashCountTask<MAX_K, DMAX_C, SuperKStorageReader>>(
                path, config, superk_storage, pinfo, i, KmDir::get().m_fof.get_i(opt->id),
                hw.get_window_size_bits(), config._kmerSize, opt->c_ab_min, opt->lz4,
                get_hist_clone(hist), opt->clear));
        }
        else
        {
          spdlog::debug("[push] - HashVecCountTask - S={}, P={}", opt->id, i);
          pool.add_task(std::make_shared<HashVecCountTask<MAX_K, DMAX_C, SuperKStorageReader>>(
            path, config, superk_storage, pinfo, i, KmDir::get().m_fof.get_i(opt->id),
            hw.get_window_size_bits(), config._kmerSize, opt->c_ab_min, opt->lz4, get_hist_clone(hist), opt->clear));
        }
      }
    }
    pool.join_all();

    if (opt->hist)
    {
      hist->merge_clones();
      HistWriter(KmDir::get().get_hist_path(opt->id), *hist, false);
    }
  }
};

template<size_t MAX_K>
struct main_merge
{
  void operator()(km_options_t options)
  {
    spdlog::info("Run with {} implementation", Kmer<MAX_K>::name());
    merge_options_t opt = std::static_pointer_cast<struct merge_options>(options);
    spdlog::debug(opt->display());
    KmDir::get().init(opt->dir, "", false);
    opt->init_vector();

    Storage* config_storage = StorageFactory(STORAGE_FILE).load(KmDir::get().m_config_storage);
    LOCAL(config_storage);
    Configuration config = Configuration();
    config.load(config_storage->getGroup("gatb"));

    if (opt->m_ab_float)
    {
      std::vector<hist_t> hist;
      for (auto id : KmDir::get().m_fof)
      {
        hist.push_back(HistReader<8192>(KmDir::get().get_hist_path(std::get<0>(id))).get());
      }
      opt->m_ab_min_vec = compute_merge_thresholds(hist, opt->m_ab_min_f, KmDir::get().get_merge_th_path());
    }

    HashWindow hw(KmDir::get().m_hash_win);

    TaskPool pool(opt->nb_threads);

    std::vector<uint32_t> ab_vec(KmDir::get().m_fof.size(), opt->m_ab_min);
    for (size_t i=0; i<config._nb_partitions; i++)
    {
      if (opt->partition_id != -1)
        if (static_cast<size_t>(opt->partition_id) != i)
          continue;
      if (opt->count_format == COUNT_FORMAT::KMER)
      {
        spdlog::debug("[push] - KmerMergeTask - P={}", i);
        pool.add_task(std::make_shared<KmerMergeTask<MAX_K, DMAX_C>>(
          i, ab_vec, config._kmerSize, opt->r_min, opt->save_if, opt->lz4, opt->mode, opt->format));
      }
      else
      {
        spdlog::debug("[push] - HashMergeTask - P={}", i);
        pool.add_task(std::make_shared<HashMergeTask<DMAX_C>>(
          i, ab_vec, opt->r_min, opt->save_if, opt->lz4, opt->mode, opt->format, hw, false, 0));
      }
    }
    pool.join_all();
  }
};

template<size_t MAX_K>
struct main_format
{
  void operator()(km_options_t options)
  {
    spdlog::info("Run with {} implementation", Kmer<MAX_K>::name());
    format_options_t opt = std::static_pointer_cast<struct format_options>(options);
    spdlog::debug(opt->display());
    KmDir::get().init(opt->dir, "", false);

    HashWindow hw(KmDir::get().m_hash_win);

    Storage* config_storage = StorageFactory(STORAGE_FILE).load(KmDir::get().m_config_storage);
    LOCAL(config_storage);
    Configuration config = Configuration();
    config.load(config_storage->getGroup("gatb"));

    if (opt->from_hash)
    {
      std::vector<vmr_t<8192>> files;
      std::vector<int> fds;
      std::vector<std::mutex> mutex(config._nb_partitions);
      for (size_t p=0; p<config._nb_partitions; p++)
      {
        fds.push_back(open(KmDir::get().get_matrix_path(p, MODE::BFT, FORMAT::BIN, COUNT_FORMAT::HASH, false).c_str(), O_RDONLY));
      }
      TaskPool pool(opt->nb_threads);

      if (opt->id == "all")
      {
        for (auto& id : KmDir::get().m_fof)
        {
          std::string sid = std::get<0>(id);
          uint32_t file_id = KmDir::get().m_fof.get_i(sid);
          pool.add_task(std::make_shared<FormatTask>(
            fds, mutex, opt->out_format, hw.bloom_size(), file_id, config._nb_partitions,
            config._kmerSize, opt->clear));
        }
      }
      pool.join_all();
    }
    else if (opt->from_vec)
    {
      TaskPool pool(opt->nb_threads);
      for (auto& id : KmDir::get().m_fof)
      {
        spdlog::debug("[push] - FormatVectorTask - S={}", std::get<0>(id));
        pool.add_task(std::make_shared<FormatVectorTask>(
          std::get<0>(id), opt->out_format, hw.bloom_size(), config._nb_partitions, opt->lz4,
          config._kmerSize, opt->clear));
      }
      pool.join_all();
    }
  }
};

template<size_t MAX_K>
struct main_dump
{
  void operator()(km_options_t options)
  {
    spdlog::info("Run with {} implementation", Kmer<MAX_K>::name());
    dump_options_t opt = std::static_pointer_cast<struct dump_options>(options);
    spdlog::debug(opt->display());

    KM_FILE km_file = get_km_file_type(opt->input);

    if (km_file == KM_FILE::KMER)
    {
      KmerReader kr(opt->input);
      if (opt->output == "stdout")
        kr.template write_as_text<MAX_K, DMAX_C>(std::cout);
      else
      {
        std::ofstream out(opt->output); check_fstream_good(opt->output, out);
        kr.template write_as_text<MAX_K, DMAX_C>(out);
      }
    }
    else if (km_file == KM_FILE::HASH)
    {
      HashReader<DMAX_C, 32768> hr(opt->input);
      if (opt->output == "stdout")
        hr.write_as_text(std::cout);
      else
      {
        std::ofstream out(opt->output); check_fstream_good(opt->output, out);
        hr.write_as_text(out);
      }
    }
    else if (km_file == KM_FILE::MATRIX)
    {
      MatrixReader mr(opt->input);
      if (opt->output == "stdout")
        mr.template write_as_text<MAX_K, DMAX_C>(std::cout);
      else
      {
        std::ofstream out(opt->output); check_fstream_good(opt->output, out);
        mr.template write_as_text<MAX_K, DMAX_C>(out);
      }
    }
    else if (km_file == KM_FILE::MATRIX_HASH)
    {
      MatrixHashReader mhr(opt->input);
      if (opt->output == "stdout")
        mhr.template write_as_text<DMAX_C>(std::cout);
      else
      {
        std::ofstream out(opt->output); check_fstream_good(opt->output, out);
        mhr.template write_as_text<DMAX_C>(out);
      }
    }
    else if (km_file == KM_FILE::PAMATRIX)
    {
      PAMatrixReader pr(opt->input);
      if (opt->output == "stdout")
        pr.template write_as_text<MAX_K>(std::cout);
      else
      {
        std::ofstream out(opt->output); check_fstream_good(opt->output, out);
        pr.template write_as_text<MAX_K>(out);
      }
    }
    else if (km_file == KM_FILE::PAMATRIX_HASH)
    {
      PAHashMatrixReader phr(opt->input);
      if (opt->output == "stdout")
        phr.write_as_text(std::cout);
      else
      {
        std::ofstream out(opt->output); check_fstream_good(opt->output, out);
        phr.write_as_text(out);
      }
    }
    else if (km_file == KM_FILE::HIST)
    {
      HistReader hr(opt->input);
      if (opt->output == "stdout")
        hr.write_as_text(std::cout, false);
      else
      {
        std::ofstream out(opt->output); check_fstream_good(opt->output, out);
        hr.write_as_text(out, false);
      }
    }
    else
    {
      throw IOError(fmt::format("KM_FILE::{} doesn't support text conversion.",
                                km_file_to_str(km_file)));
    }
  }
};

template<size_t MAX_K>
struct main_combine
{
  void operator()(km_options_t options)
  {
    spdlog::info("Run with {} implementation", Kmer<MAX_K>::name());
    combine_options_t opt = std::static_pointer_cast<struct combine_options>(options);
    spdlog::debug(opt->display());

    auto parse_mode = [](const std::string& p) {
      std::ifstream inf(fmt::format("{}/options.txt", p), std::ios::in);
      std::string line; std::getline(inf, line);

      MODE m; COUNT_FORMAT c;
      auto v = bc::utils::split(line, ',');
      for (auto& e : v)
      {
        auto vv = bc::utils::split(e, '=');
        auto entry = bc::utils::trim(vv[0]);

        if (entry == "mode")
          m = str_to_mode(bc::utils::trim(vv[1]));
        else if (entry == "count_format")
          c = str_to_cformat(bc::utils::trim(vv[1]));
      }

      return std::make_tuple(m, c);
    };

    Timer timer;

    auto [m, c] = parse_mode(opt->runs[0]);

    for (auto& cc : opt->runs)
      spdlog::info(cc);
    TaskPool pool(opt->nb_threads);

    if (m == MODE::COUNT && c == COUNT_FORMAT::KMER)
    {
      MatrixMerger<MAX_K, DMAX_C> mm(opt->runs, opt->output, opt->cpr);
      mm.exec(pool);

    }
    else if (m == MODE::PA && c == COUNT_FORMAT::KMER)
    {
      MatrixMerger<MAX_K, 1> mm(opt->runs, opt->output, opt->cpr);
      mm.exec(pool);
    }
    else if (m == MODE::COUNT && c == COUNT_FORMAT::HASH)
    {
      MatrixMerger<1, DMAX_C> mm(opt->runs, opt->output, opt->cpr);
      mm.exec(pool);
    }
    else if (m == MODE::PA && c == COUNT_FORMAT::HASH)
    {
      MatrixMerger<1, 1> mm(opt->runs, opt->output, opt->cpr);
      mm.exec(pool);
    }
    else
    {
      spdlog::debug("mode = {}, count = {}", mode_to_str(m), cformat_to_str(c));
      throw InputError(
        fmt::format("{}: matrix format not supported by 'kmtricks combine'.", opt->runs[0]));
    }

    spdlog::info("Done in {}. New matrix is located at {}.", timer.formatted(),  opt->output);
  }
};


template<size_t MAX_K>
struct main_agg
{
  void operator()(km_options_t options)
  {
    spdlog::info("Run with {} implementation", Kmer<MAX_K>::name());
    agg_options_t opt = std::static_pointer_cast<struct agg_options>(options);
    spdlog::debug(opt->display());

    KmDir::get().init(opt->dir, "", false);
    Storage* config_storage = StorageFactory(STORAGE_FILE).load(KmDir::get().m_config_storage);
    LOCAL(config_storage);
    Configuration config = Configuration();
    config.load(config_storage->getGroup("gatb"));

    auto check_paths = [](const std::vector<std::string>& paths) {
      std::vector<std::string> ret;
      for (auto& p : paths)
      {
        if (fs::exists(p))
          ret.push_back(p);
      }

      if (ret.empty())
        throw IOError("No files found for these parameters.");
      return ret;
    };

    if (opt->count == "kmer")
    {
      std::vector<std::string> paths = KmDir::get().get_count_part_paths(
        opt->id, config._nb_partitions, opt->lz4_in, KM_FILE::KMER);
      paths = check_paths(paths);

      if (opt->sorted)
      {
        KmerFileMerger<MAX_K, DMAX_C> kfm(paths, config._kmerSize);
        if (opt->format == "text")
        {
          if (opt->no_count)
            opt->output == "stdout" ? kfm.write_kmers(std::cout) : kfm.write_kmers(opt->output);
          else
            opt->output == "stdout" ? kfm.write_as_text(std::cout) : kfm.write_as_text(opt->output);
        }
        else
          kfm.write_as_bin(opt->output, opt->lz4);
      }
      else
      {
        KmerFileAggregator<MAX_K, DMAX_C> kfa(paths, config._kmerSize);
        if (opt->format == "text")
        {
          if (opt->no_count)
            opt->output == "stdout" ? kfa.write_kmers(std::cout) : kfa.write_kmers(opt->output);
          else
            opt->output == "stdout" ? kfa.write_as_text(std::cout) : kfa.write_as_text(opt->output);
        }
        else
          kfa.write_as_bin(opt->output, opt->lz4);
      }
    }
    else if (opt->count == "hash")
    {
      std::vector<std::string> paths = KmDir::get().get_count_part_paths(
        opt->id, config._nb_partitions, opt->lz4_in, KM_FILE::HASH);
      paths = check_paths(paths);

      HashFileAggregator<DMAX_C> hfa(paths);
      if (opt->format == "text")
        opt->output == "stdout" ? hfa.write_as_text(std::cout) : hfa.write_as_text(opt->output);
      else
        hfa.write_as_bin(opt->output, opt->lz4);
    }
    else if (opt->matrix == "kmer")
    {
      std::vector<std::string> paths = KmDir::get().get_matrix_paths(config._nb_partitions,
                                                                    MODE::COUNT, FORMAT::BIN,
                                                                    COUNT_FORMAT::KMER, opt->lz4_in);
      paths = check_paths(paths);

      if (opt->sorted)
      {
        MatrixFileMerger<MAX_K, DMAX_C> mfm(paths, config._kmerSize);
        if (opt->format == "text")
        {
          if (opt->no_count)
            opt->output == "stdout" ? mfm.write_kmers(std::cout) : mfm.write_kmers(opt->output);
          else
            opt->output == "stdout" ? mfm.write_as_text(std::cout) : mfm.write_as_text(opt->output);
        }
        else
          mfm.write_as_bin(opt->output, opt->lz4);
      }
      else
      {
        MatrixFileAggregator<MAX_K, DMAX_C> mfa(paths, config._kmerSize);
        if (opt->format == "text")
        {
          if (opt->no_count)
            opt->output == "stdout" ? mfa.write_kmers(std::cout) : mfa.write_kmers(opt->output);
          else
            opt->output == "stdout" ? mfa.write_as_text(std::cout) : mfa.write_as_text(opt->output);
        }
        else
          mfa.write_as_bin(opt->output, opt->lz4);
      }
    }
    else if (opt->matrix == "hash")
    {
      std::vector<std::string> paths = KmDir::get().get_matrix_paths(config._nb_partitions,
                                                                    MODE::COUNT, FORMAT::BIN,
                                                                    COUNT_FORMAT::HASH, opt->lz4_in);
      paths = check_paths(paths);

      MatrixHashFileAggregator<DMAX_C> mhfa(paths);
      if (opt->format == "text")
        opt->output == "stdout" ? mhfa.write_as_text(std::cout) : mhfa.write_as_text(opt->output);
      else
        mhfa.write_as_bin(opt->output, opt->lz4);
    }
    else if (opt->pa_matrix == "kmer")
    {
      std::vector<std::string> paths = KmDir::get().get_matrix_paths(config._nb_partitions,
                                                                    MODE::PA, FORMAT::BIN,
                                                                    COUNT_FORMAT::KMER, opt->lz4_in);
      paths = check_paths(paths);
      if (opt->sorted)
      {
        PAMatrixFileMerger<MAX_K> pmfm(paths, config._kmerSize);
        if (opt->format == "text")
        {
          if (opt->no_count)
            opt->output == "stdout" ? pmfm.write_as_text(std::cout) : pmfm.write_as_text(opt->output);
          else
            opt->output == "stdout" ? pmfm.write_as_text(std::cout) : pmfm.write_as_text(opt->output);
        }
        else
          pmfm.write_as_bin(opt->output, opt->lz4);
      }
      else
      {
        PAMatrixFileAggregator<MAX_K> pmfa(paths, config._kmerSize);
        if (opt->format == "text")
        {
          if (opt->no_count)
            opt->output == "stdout" ? pmfa.write_kmers(std::cout) : pmfa.write_kmers(opt->output);
          else
            opt->output == "stdout" ? pmfa.write_as_text(std::cout) : pmfa.write_as_text(opt->output);
        }
        else
          pmfa.write_as_bin(opt->output, opt->lz4);
      }
    }
    else if (opt->pa_matrix == "hash")
    {
      std::vector<std::string> paths = KmDir::get().get_matrix_paths(config._nb_partitions,
                                                                    MODE::PA, FORMAT::BIN,
                                                                    COUNT_FORMAT::HASH, opt->lz4_in);
      paths = check_paths(paths);
      PAHashMatrixFileAggregator phmfa(paths);
      if (opt->format == "text")
        opt->output == "stdout" ? phmfa.write_as_text(std::cout) : phmfa.write_as_text(opt->output);
      else
        phmfa.write_as_bin(opt->output, opt->lz4);
    }
  }
};

template<size_t MAX_K>
struct main_filter
{
  void operator()(km_options_t options)
  {
    spdlog::info("Run with {} implementation", Kmer<MAX_K>::name());
    filter_options_t opt = std::static_pointer_cast<struct filter_options>(options);

    KmDir::get().init(opt->dir, "", false);
    std::string in_config = fmt::format("{}_gatb", KmDir::get().m_config_storage);
    std::string in_repart = fmt::format("{}_gatb", KmDir::get().m_repart_storage);

    Storage* config_storage = StorageFactory(STORAGE_FILE).load(KmDir::get().m_config_storage);
    LOCAL(config_storage);
    Configuration config = Configuration();
    config.load(config_storage->getGroup("gatb"));

    std::vector<std::string> in_matrices;
    std::vector<uint32_t> partitions;

    MODE mode = MODE::COUNT;

    for (std::uint32_t p = 0; p < config._nb_partitions; p++)
    {
      std::string mp = KmDir::get().get_matrix_path(p, MODE::PA, FORMAT::BIN, COUNT_FORMAT::KMER, opt->cpr_in);
      std::string mc = KmDir::get().get_matrix_path(p, MODE::COUNT, FORMAT::BIN, COUNT_FORMAT::KMER, opt->cpr_in);

      if (fs::exists(mp))
      {
        mode = MODE::PA;
        in_matrices.push_back(mp);
        partitions.push_back(p);
      }
      else if (fs::exists(mc))
      {
        mode = MODE::COUNT;
        in_matrices.push_back(mc);
        partitions.push_back(p);
      }
    }

    if (in_matrices.empty())
      throw IOError("No files found for these parameters");

    KmDir::get().init(opt->output, opt->key, true);

    if (KmDir::get().m_fof.size() > 1)
      throw InputError("Filtering with many samples is not yet implemented. Fof must contain only one sample.");

    fs::copy(in_config, fmt::format("{}_gatb", KmDir::get().m_config_storage));
    fs::copy(in_repart, fmt::format("{}_gatb", KmDir::get().m_repart_storage));

    std::string sid = KmDir::get().m_fof.get_id(0);

    spdlog::info("Key = {}", sid);
    spdlog::info("Compute super-k-mers (process {} partition(s))...", partitions.size());
    SuperKTask<MAX_K> superk_task(sid, true, partitions);
    superk_task.exec();

    sk_storage_t superk_storage = std::make_shared<SuperKStorageReader>(
      KmDir::get().get_superk_path(sid));
    parti_info_t pinfo = std::make_shared<PartiInfo<5>>(KmDir::get().get_superk_path(sid));

    TaskPool pool(opt->nb_threads);

    std::size_t amin = std::get<2>(*(KmDir::get().m_fof.begin()));
    amin = (amin == 0) ? opt->c_ab_min : amin;

    spdlog::info("Count partitions...");
    for (auto&& i : partitions)
    {
      KmDir::get().init_one_part(i);
      std::string p = KmDir::get().get_count_part_path(sid, i, true, KM_FILE::KMER);
      uint32_t id = KmDir::get().m_fof.get_i(sid);

      pool.add_task(std::make_shared<CountTask<MAX_K, DMAX_C, SuperKStorageReader>>(
        p, config, superk_storage, pinfo, i, id, config._kmerSize, amin, true, nullptr, false
      ));
    }
    pool.join_all();

    std::vector<std::string> out_matrices;
    std::vector<std::string> in_kmers;
    std::vector<std::string> out_kmers;
    std::vector<std::string> vecs;

    for (auto&& p : partitions)
    {
      out_matrices.push_back(
        KmDir::get().get_matrix_path(p, mode, FORMAT::BIN, COUNT_FORMAT::KMER, opt->cpr_out));
      in_kmers.push_back(
        KmDir::get().get_count_part_path(sid, p, true, KM_FILE::KMER));
      out_kmers.push_back(
        KmDir::get().get_count_part_path(fmt::format("{}_absent", sid), p, opt->cpr_out, KM_FILE::KMER));
      vecs.push_back(
        fmt::format("{}/{}.vec", KmDir::get().m_matrix_storage, p));
    }

    std::tuple<bool, bool, bool> out_types = std::make_tuple(opt->with_vector, opt->with_matrix, opt->with_kmer);

    spdlog::info("Filtering...");
    MatrixFilter<MAX_K, DMAX_C> mf(in_matrices, in_kmers, out_matrices, out_kmers, vecs, opt->cpr_out, mode == MODE::COUNT, opt->nb_threads, out_types);
    mf.exec();

    for (std::size_t i = 0; i < partitions.size(); ++i)
    {
      if (opt->with_kmer)
      {
        fs::rename(out_kmers[i], KmDir::get().get_count_part_path(sid, partitions[i], opt->cpr_out, KM_FILE::KMER));
      }
      else
      {
        fs::remove(KmDir::get().get_count_part_path(sid, partitions[i], true, KM_FILE::KMER));
      }
    }
  }
};

#ifdef WITH_HOWDE
template<size_t MAX_K>
struct main_index
{
  void operator()(km_options_t options)
  {
    spdlog::info("Run with {} implementation", Kmer<MAX_K>::name());
    index_options_t opt = std::static_pointer_cast<struct index_options>(options);
    spdlog::debug(opt->display());
    KmDir::get().init(opt->dir, "", false);
    spdlog::info("Compute tree topology...");
    std::string bf_list = KmDir::get().get_bf_list_path();
    {
      std::ofstream out(bf_list, std::ios::out); check_fstream_good(bf_list, out);
      for (auto id: KmDir::get().m_fof)
        out << fs::absolute(fs::path(
          KmDir::get().get_filter_path(std::get<0>(id), OUT_FORMAT::HOWDE))).string() << "\n";
    }

    std::string index = KmDir::get().get_index_path();
    std::string howde_index_str = fmt::format("cluster --list={} --out={}",
                                              bf_list, index);

    if (opt->upper != 0)
      howde_index_str += fmt::format(" {}..{}", opt->lower, opt->upper);
    else
      howde_index_str += fmt::format(" --bits={}", opt->bits);

    if (opt->cull > 0)
      howde_index_str += fmt::format(" --cull={}", opt->cull);

    if (opt->cull2)
      howde_index_str += " --cull";

    if (opt->cullsd > 0)
      howde_index_str += fmt::format(" --cull={}sd", opt->cullsd);
    std::vector<std::string> howde_index = bc::utils::split(howde_index_str, ' ');

    char** arr = new char*[howde_index.size()+1];
    arr[howde_index.size()] = nullptr;
    for (size_t i=0; i<howde_index.size(); i++)
      arr[i] = strdup(howde_index.at(i).c_str());

    ClusterCommand cluster_cmd("cluster");
    cluster_cmd.parse(howde_index.size(), arr);
    cluster_cmd.execute();

    for (size_t i=0; i<howde_index.size(); i++)
      free(arr[i]);

    delete[] arr;

    spdlog::info("Build index...");
    std::stringstream ss;
    ss << "build ";
    ss << KmDir::get().get_index_path() << " ";
    if (opt->howde) ss << "--howde ";
    if (opt->allsome) ss << "--allsome ";
    if (opt->determined) ss << "--determined ";
    if (opt->brief) ss << "--determined,brief ";
    if (opt->uncompressed) ss << "--uncompressed ";
    if (opt->rrr) ss << "--rrr ";
    if (opt->roar) ss << "--roar ";
    std::string howde_build_str = ss.str();
    std::vector<std::string> howde_build = bc::utils::split(howde_build_str, ' ');

    char** arr2 = new char*[howde_build.size()+1];
    arr2[howde_build.size()] = nullptr;
    for (size_t i=0; i<howde_build.size(); i++)
      arr2[i] = strdup(howde_build.at(i).c_str());

    BuildSBTCommand build_cmd("build");
    build_cmd.parse(howde_build.size(), arr2);
    auto path = fs::current_path();

    fs::current_path(KmDir::get().m_index_storage);
    build_cmd.execute();

    fs::current_path(path);
    for (size_t i=0; i<howde_build.size(); i++)
      free(arr2[i]);
    delete[] arr2;
  }
};

template<size_t MAX_K>
struct main_query
{
  void operator()(km_options_t options)
  {
    spdlog::info("Run with {} implementation", Kmer<MAX_K>::name());
    query_options_t opt = std::static_pointer_cast<struct query_options>(options);
    spdlog::debug(opt->display());

    KmDir::get().init(opt->dir, "", false);
    Storage* config_storage = StorageFactory(STORAGE_FILE).load(KmDir::get().m_config_storage);
    LOCAL(config_storage);
    Configuration config = Configuration();
    config.load(config_storage->getGroup("gatb"));

    std::string index_path;
    for (auto& p : fs::directory_iterator(KmDir::get().m_index_storage))
    {
      if (p.path().string().find(".sbt") != std::string::npos)
      {
        index_path = p.path().string();
        break;
      }
    }

    if (index_path.empty())
      throw IOError("Index not found.");

    if (opt->output != "stdout")
      opt->output = fs::absolute(fs::path(opt->output));

    opt->query = fs::absolute(fs::path(opt->query));

    std::stringstream ss;
    ss << "queryKm ";
    ss << "--tree=" << index_path << " ";
    ss << opt->query << " ";
    ss << "--repart=" << fmt::format("{}_gatb/repartition.minimRepart", KmDir::get().m_repart_storage) << " ";
    ss << "--win=" << KmDir::get().m_hash_win << " ";
    ss << "--z=" << opt->z << " ";
    ss << "--threshold=" << opt->threshold << " ";
    ss << "--threshold-shared-positions=" << opt->threshold_shared_positions << " ";
    if (opt->check) ss << "--consistencycheck ";
    if (opt->nodetail) ss << "--no-detail ";
    if (opt->output != "stdout") ss << "--out=" << opt->output;

    std::string howde_query_str = ss.str();
    std::vector<std::string> howde_query = bc::utils::split(howde_query_str, ' ');

    char** arr = new char*[howde_query.size()+1];
    arr[howde_query.size()] = nullptr;
    for (size_t i=0; i<howde_query.size(); i++)
      arr[i] = strdup(howde_query.at(i).c_str());

    QueryCommand query_cmd("queryKm");
    query_cmd.parse(howde_query.size(), arr);
    auto path = fs::current_path();

    fs::current_path(KmDir::get().m_index_storage);
    query_cmd.execute();

    fs::current_path(path);
    for (size_t i=0; i<howde_query.size(); i++)
      free(arr[i]);
    delete[] arr;
  }
};

#endif
};
