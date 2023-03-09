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
#include <cmath>
#include <functional>
#include <filesystem>

#include <gatb/gatb_core.hpp>
#include <gatb/kmer/impl/RepartitionAlgorithm.hpp>
#include <gatb/kmer/impl/ConfigurationAlgorithm.hpp>
#include <gatb/kmer/impl/SortingCountAlgorithm.cpp>

#include <kmtricks/io/fof.hpp>
#include <kmtricks/kmdir.hpp>
#include <kmtricks/gatb/count_processor.hpp>
#include <kmtricks/gatb/sorting_count.hpp>
#include <kmtricks/gatb/fill_partitions.hpp>
#include <kmtricks/merge.hpp>
#include <kmtricks/hash.hpp>
#include <kmtricks/howde_utils.hpp>
#include <kmtricks/gatb/gatb_utils.hpp>
#include <kmtricks/itask.hpp>
#include <kmtricks/repartition.hpp>

#ifdef WITH_PLUGIN
#include <kmtricks/plugin_manager.hpp>
#endif

namespace km {

namespace fs = std::filesystem;
using parti_info_t = std::shared_ptr<PartiInfo<5>>;

template<size_t span>
class ConfigTask : public ITask
{
public:
  ConfigTask(const std::string& path, IProperties* props, uint64_t bloom_size,
             uint32_t partitions)
    : ITask(0), m_path(path), m_props(props), m_bloom_size(bloom_size), m_nb_partitions(partitions)
  {}

  void preprocess() {}
  void postprocess() {}
  void exec()
  {
    spdlog::debug("[exec] - ConfigTask");
    spdlog::info("{} samples found ({} read files).", KmDir::get().m_fof.size(), KmDir::get().m_fof.total());
    IBank* bank = Bank::open(KmDir::get().m_fof.get_all()); LOCAL(bank);
    Storage* config_storage =
      StorageFactory(STORAGE_FILE).create(KmDir::get().m_config_storage, true, false);
    LOCAL(config_storage);

    ConfigurationAlgorithm<span> config_alg(bank, m_props);
    config_alg.execute();
    Configuration config = config_alg.getConfiguration();

    if (m_nb_partitions != 0)
      config._nb_partitions = m_nb_partitions;
    if (config._nb_partitions < 4)
      config._nb_partitions = 4;

    spdlog::info("Use {} partitions.", config._nb_partitions);

    config.save(config_storage->getGroup("gatb"));
    HashWindow hw(m_bloom_size, config._nb_partitions, config._minim_size);
    hw.serialize(KmDir::get().m_hash_win);

    spdlog::debug("[done] - ConfigTask");
  }

private:
  std::string m_path;
  IProperties* m_props;
  uint64_t m_bloom_size;
  uint32_t m_nb_partitions;
};

void check_repart_compatibility(Configuration& c1, Configuration& c2,
                                const std::string& d1)
{
  static std::string err_temp = "Unable to use repartition from {}, {} differ.";

  if (c1._kmerSize != c2._kmerSize)
    throw PipelineError(fmt::format(err_temp, d1, "kmer sizes"));
  if (c1._minim_size != c2._minim_size)
    throw PipelineError(fmt::format(err_temp, d1, "minimizer sizes"));
  if (c1._nb_partitions != c2._nb_partitions)
    throw PipelineError(fmt::format(err_temp, d1, "numbers of partitions"));
}

template<size_t span>
class RepartTask : public ITask
{
public:
  RepartTask(const std::string& path, const std::string& from = "")
    : ITask(1), m_path(path), m_from(from) {}

  void preprocess() {}
  void postprocess()
  {
    if (m_minim_size <= 12)
    {
      std::vector<std::string> paths = KmDir::get().get_minim_paths(m_nb_parts);
      Repartition repart(fmt::format("{}_gatb/repartition.minimRepart", KmDir::get().m_repart_storage));
      repart.write_minimizers(paths, m_minim_size);
    }
  }

  void exec()
  {
    spdlog::debug("[exec] - RepartTask");

    Storage* config_storage =
      StorageFactory(STORAGE_FILE).load(KmDir::get().m_config_storage);
    LOCAL(config_storage);

    Configuration config = Configuration();
    config.load(config_storage->getGroup("gatb"));
    m_nb_parts = config._nb_partitions;
    m_minim_size = config._minim_size;

    if (m_from.empty())
    {
      Fof fof(m_path);
      IBank* bank = Bank::open(fof.get_all()); LOCAL(bank);
      Storage* rep_store =
        StorageFactory(STORAGE_FILE).create(KmDir::get().m_repart_storage, true, false);

      LOCAL(rep_store);

      RepartitorAlgorithm<span> repartition(
        bank, rep_store->getGroup("repartition"), config, 1);
      repartition.execute();
    }
    else
    {
      Storage* fc_store = StorageFactory(STORAGE_FILE).load(
        fmt::format("{}/{}", m_from, "config"));
      LOCAL(fc_store);
      Configuration fc_config;
      fc_config.load(fc_store->getGroup("gatb"));

      check_repart_compatibility(config, fc_config, m_from);

      fs::copy(
        fmt::format("{}/{}", m_from, "repartition_gatb"),
        fmt::format("{}/{}", KmDir::get().m_root, "repartition_gatb"),
        fs::copy_options::recursive | fs::copy_options::overwrite_existing);
    }

    spdlog::debug("[done] - RepartTask");
  }

private:
  std::string m_path;
  std::string m_from;
  int m_cores;
  uint32_t m_nb_parts {0};
  uint32_t m_minim_size {0};
};

template<size_t span>
class SuperKTask : public ITask
{
public:
  SuperKTask(const std::string& sample_id, bool lz4, std::vector<uint32_t>& partitions)
    : ITask(2), m_sample_id(sample_id), m_lz4(lz4), m_partitions(partitions) {}

  void preprocess() {}

  void postprocess()
  {
    this->exec_callback();
    this->m_finish = true;
    this->m_running = false;
  }

  void exec()
  {
    spdlog::debug("[exec] - SuperKTask - S={}", m_sample_id);
    this->m_running = true;

    IBank* bank = Bank::open(KmDir::get().m_fof.get_files(m_sample_id)); LOCAL(bank);
    Storage* config_storage = StorageFactory(STORAGE_FILE).load(KmDir::get().m_config_storage);
    Storage* repart_storage = StorageFactory(STORAGE_FILE).load(KmDir::get().m_repart_storage);
    LOCAL(config_storage); LOCAL(repart_storage);

    Configuration config = Configuration();
    config.load(config_storage->getGroup("gatb"));

    Repartitor repartitor(repart_storage->getGroup("repartition"));

    std::unordered_set<int> pset;
    for (auto& p : m_partitions)
    {
      pset.insert(p);
    }
    SuperKStorageWriter* superk_storage = new SuperKStorageWriter(
      KmDir::get().get_superk_path(m_sample_id), "skp", config._nb_partitions, m_lz4, pset);

    typedef typename ::Kmer<span>::ModelCanonical ModelCanonical;
    typedef typename ::Kmer<span>::template ModelMinimizer <ModelCanonical> Model;

    uint32_t* freq_order = nullptr;
    Model model(config._kmerSize, config._minim_size,
                typename ::Kmer<span>::ComparatorMinimizerFrequencyOrLex(), freq_order);

    Iterator<Sequence>* itSeq = bank->iterator(); LOCAL(itSeq);
    BankStats bank_stats;
    PartiInfo<5> pinfo (config._nb_partitions, config._minim_size);

    IteratorListener* progress(new ProgressSynchro(
                               new IteratorListener(),
                               System::thread().newSynchronizer()));
    LOCAL(progress);
    progress->init();
    {
      auto fill_partitions = KmFillPartitions<span>(model,
                                                        1,
                                                        0,
                                                        config._nb_partitions,
                                                        config._nb_cached_items_per_core_per_part,
                                                        progress,
                                                        bank_stats,
                                                        nullptr,
                                                        repartitor,
                                                        pinfo,
                                                        superk_storage);

      for (itSeq->first(); !itSeq->isDone(); itSeq->next())
      {
        fill_partitions(itSeq->item());
      }
      itSeq->finalize();
    }

    progress->finish();
    superk_storage->SaveInfoFile(KmDir::get().get_superk_path(m_sample_id));
    delete superk_storage;
    pinfo.saveInfoFile(KmDir::get().get_superk_path(m_sample_id));
    dump_pinfo(&pinfo, config._nb_partitions, KmDir::get().get_pinfos_path(m_sample_id));
    spdlog::debug("[done] - SuperKTask - S={}", m_sample_id);
  }

private:
  std::string m_sample_id;
  bool m_lz4;
  std::vector<uint32_t>& m_partitions;
};

template<size_t span, size_t MAX_C, typename Storage>
class CountTask : public ITask
{
  using storage_t = std::shared_ptr<Storage>;
public:
  CountTask(const std::string& path,
            Configuration& config,
            storage_t superk_storage,
            parti_info_t pinfo,
            uint32_t part_id, uint32_t sample_id,
            uint32_t kmer_size, uint32_t abundance_min, bool lz4,
            hist_t hist = nullptr, bool clear = false)
    : ITask(3, clear),
      m_path(path),
      m_config(config),
      m_superk_storage(superk_storage),
      m_pinfo(pinfo),
      m_part_id(part_id),
      m_sample_id(sample_id),
      m_kmer_size(kmer_size),
      m_ab_min(abundance_min),
      m_lz4(lz4),
      m_hist(hist)
   {
   }

  void preprocess() {}
  void postprocess()
  {
    if (this->m_clear)
    {
      m_superk_storage->closeFile(m_part_id);
      //m_superk_storage->eraseFile(m_part_id);
      Eraser::get().erase(m_superk_storage->getFileName(m_part_id));
    }
    this->m_finish = true;
    this->exec_callback();
  }

  void exec()
  {
    spdlog::debug("[exec] - CountTask - S={}, P={}", KmDir::get().m_fof.get_id(m_sample_id), m_part_id);

    MemAllocator pool(1);
    pool.reserve(get_required_memory<span>(m_pinfo->getNbKmer(m_part_id)));
    kw_t<8192> writer = std::make_shared<KmerWriter<8192>>(m_path,
                                                           m_kmer_size,
                                                           requiredC<MAX_C>::value/8,
                                                           m_sample_id,
                                                           m_part_id,
                                                           m_lz4);

    KmerCountProcessor<span, MAX_C>* processor(new KmerCountProcessor<span, MAX_C>(m_kmer_size,
                                                                                    m_ab_min,
                                                                                    writer,
                                                                                    m_hist));

    KmerPartCounter<Storage, span> partition_counter(processor, m_pinfo.get(), m_part_id,
                                                     m_kmer_size, pool, m_superk_storage.get());

    partition_counter.execute();
    pool.free_all();
    delete processor;
    spdlog::debug("[done] - CountTask - S={}, P={}", KmDir::get().m_fof.get_id(m_sample_id), m_part_id);
  }
private:
  std::string m_path;
  Configuration& m_config;
  //Storage& m_superk_storage;
  //PartiInfo<5>& m_pinfo;
  storage_t m_superk_storage;
  parti_info_t m_pinfo;
  uint32_t m_part_id;
  uint32_t m_sample_id;
  uint32_t m_kmer_size;
  uint32_t m_ab_min;
  bool m_lz4;
  hist_t m_hist;
};

template<size_t span, size_t MAX_C, typename Storage>
class HashCountTask : public ITask
{
  using storage_t = std::shared_ptr<Storage>;
public:
  HashCountTask(const std::string& path,
                Configuration& config,
                storage_t superk_storage,
                parti_info_t pinfo,
                uint32_t part_id, uint32_t sample_id, uint64_t window,
                uint32_t kmer_size, uint32_t abundance_min, bool lz4,
                hist_t hist = nullptr, bool clear = false)
    : ITask(3, clear),
      m_path(path),
      m_config(config),
      m_superk_storage(superk_storage),
      m_pinfo(pinfo),
      m_part_id(part_id),
      m_sample_id(sample_id),
      m_window(window),
      m_kmer_size(kmer_size),
      m_ab_min(abundance_min),
      m_lz4(lz4),
      m_hist(hist)
   {
   }

  void preprocess() {}
  void postprocess()
  {
    if (this->m_clear)
    {
      m_superk_storage->closeFile(m_part_id);
      Eraser::get().erase(m_superk_storage->getFileName(m_part_id));
    }
    this->m_finish = true;
    this->exec_callback();
  }

  void exec()
  {
    spdlog::debug("[exec] - HashCountTask - S={}, P={}", KmDir::get().m_fof.get_id(m_sample_id), m_part_id);

    size_t nbk = m_pinfo->getNbKmer(m_part_id);

    uint64_t req_mem = get_required_memory_hash<span>(nbk);

    hw_t<MAX_C, 32768> writer = std::make_shared<HashWriter<MAX_C, 32768>>(m_path,
                                                                             requiredC<MAX_C>::value/8,
                                                                             m_sample_id,
                                                                             m_part_id,
                                                                             m_lz4);

    HashCountProcessor<span, MAX_C, 32768>* processor(new HashCountProcessor<span, MAX_C, 32768>(m_kmer_size,
                                                                                                 m_ab_min,
                                                                                                 writer,
                                                                                                 m_hist));

    if (nbk > 0)
    {
      MemAllocator pool(1);
      pool.reserve(req_mem);

      HashPartCounter<Storage, span> partition_counter(processor, m_pinfo.get(), m_part_id, m_kmer_size,
                                                       pool, m_superk_storage.get(), m_window);

      partition_counter.execute();
      pool.free_all();
    }

    delete processor;

    spdlog::debug("[done] - HashCountTask - S={}, P={}", KmDir::get().m_fof.get_id(m_sample_id), m_part_id);
  }

private:
  std::string m_path;
  Configuration& m_config;
  storage_t m_superk_storage;
  parti_info_t m_pinfo;
  uint32_t m_part_id;
  uint32_t m_sample_id;
  uint32_t m_kmer_size;
  uint64_t m_window;
  uint32_t m_ab_min;
  hist_t m_hist;
  bool m_lz4;
};

template<size_t span, size_t MAX_C, typename Storage>
class HashVecCountTask : public ITask
{
  using storage_t = std::shared_ptr<Storage>;
public:
  HashVecCountTask(const std::string& path,
                Configuration& config,
                storage_t superk_storage,
                parti_info_t pinfo,
                uint32_t part_id, uint32_t sample_id, uint64_t window,
                uint32_t kmer_size, uint32_t abundance_min, bool lz4,
                hist_t hist = nullptr, bool clear = false)
    : ITask(3, clear),
      m_path(path),
      m_config(config),
      m_superk_storage(superk_storage),
      m_pinfo(pinfo),
      m_part_id(part_id),
      m_sample_id(sample_id),
      m_window(window),
      m_kmer_size(kmer_size),
      m_ab_min(abundance_min),
      m_lz4(lz4),
      m_hist(hist)
   {
   }

  void preprocess() {}
  void postprocess()
  {
    if (this->m_clear)
    {
      m_superk_storage->closeFile(m_part_id);
      //m_superk_storage->eraseFile(m_part_id);
      Eraser::get().erase(m_superk_storage->getFileName(m_part_id));
    }
    this->m_finish = true;
    this->exec_callback();
  }

  void exec()
  {
    spdlog::debug("[exec] - HashVecCountTask - S={}, P={}",KmDir::get().m_fof.get_id(m_sample_id), m_part_id);

    size_t nbk = m_pinfo->getNbKmer(m_part_id);

    bvw_t<8192> writer = std::make_shared<BitVectorWriter<8192>>(m_path,
                                                                m_window,
                                                                0,
                                                                m_part_id,
                                                                m_lz4);

    HashVecProcessor<span>* processor(new HashVecProcessor<span>(m_kmer_size,
                                                                 m_ab_min,
                                                                 writer,
                                                                 m_hist,
                                                                 m_window));

    if (nbk > 0)
    {
      MemAllocator pool(1);
      pool.reserve(get_required_memory_hash<span>(nbk));

      HashPartCounter<Storage, span> partition_counter(processor, m_pinfo.get(), m_part_id, m_kmer_size,
                                            pool, m_superk_storage.get(), m_window);

      partition_counter.execute();
      pool.free_all();
    }

    delete processor;
    spdlog::debug("[done] - HashVecCountTask - S={}, P={}",KmDir::get().m_fof.get_id(m_sample_id), m_part_id);
  }
private:
  std::string m_path;
  Configuration& m_config;
  //Storage& m_superk_storage;
  //PartiInfo<5>& m_pinfo;
  storage_t m_superk_storage;
  parti_info_t m_pinfo;
  uint32_t m_part_id;
  uint32_t m_sample_id;
  uint32_t m_kmer_size;
  uint64_t m_window;
  uint32_t m_ab_min;
  bool m_lz4;
  hist_t m_hist;
};

template<size_t span, size_t MAX_C, typename Storage>
class KffCountTask : public ITask
{
  using storage_t = std::shared_ptr<Storage>;
public:
  KffCountTask(const std::string& path,
            Configuration& config,
            storage_t superk_storage,
            parti_info_t pinfo,
            uint32_t part_id, uint32_t sample_id,
            uint32_t kmer_size, uint32_t abundance_min,
            hist_t hist = nullptr, bool clear = false)
    : ITask(3, clear),
      m_path(path),
      m_config(config),
      m_superk_storage(superk_storage),
      m_pinfo(pinfo),
      m_part_id(part_id),
      m_sample_id(sample_id),
      m_kmer_size(kmer_size),
      m_ab_min(abundance_min),
      m_hist(hist)
   {
   }

  void preprocess() {}
  void postprocess()
  {
    if (this->m_clear)
    {
      m_superk_storage->closeFile(m_part_id);
      //m_superk_storage->eraseFile(m_part_id);
      Eraser::get().erase(m_superk_storage->getFileName(m_part_id));
    }
    this->m_finish = true;
    this->exec_callback();
  }

  void exec()
  {
    spdlog::debug("[exec] - KffCountTask - S={}, P={}", KmDir::get().m_fof.get_id(m_sample_id), m_part_id);

    MemAllocator pool(1);
    pool.reserve(get_required_memory<span>(m_pinfo->getNbKmer(m_part_id)));
    kff_w_t<DMAX_C> writer = std::make_shared<KffWriter<MAX_C>>(m_path, m_kmer_size);

    KffCountProcessor<span, DMAX_C>* processor(new KffCountProcessor<span, MAX_C>(m_kmer_size,
                                                                                  m_ab_min,
                                                                                  writer,
                                                                                  m_hist));

    KmerPartCounter<Storage, span> partition_counter(processor, m_pinfo.get(), m_part_id, m_kmer_size,
                                                     pool, m_superk_storage.get());

    partition_counter.execute();
    pool.free_all();
    delete processor;

    spdlog::debug("[done] - KffCountTask - S={}, P={}", KmDir::get().m_fof.get_id(m_sample_id), m_part_id);
  }
private:
  std::string m_path;
  Configuration& m_config;
  storage_t m_superk_storage;
  parti_info_t m_pinfo;
  uint32_t m_part_id;
  uint32_t m_sample_id;
  uint32_t m_kmer_size;
  uint32_t m_ab_min;
  hist_t m_hist;
};

template<size_t span, size_t MAX_C>
class KmerMergeTask : public ITask
{
public:
  KmerMergeTask(uint32_t partition_id,
                std::vector<uint32_t>& ab_vec,
                uint32_t kmer_size,
                uint32_t recurrence_min,
                uint32_t save_if,
                bool lz4,
                MODE mode,
                FORMAT format,
                bool clear = false)
    : ITask(4, clear), m_part_id(partition_id), m_ab_vec(ab_vec), m_kmer_size(kmer_size),
      m_rec_min(recurrence_min), m_save_if(save_if), m_lz4(lz4), m_mode(mode), m_format(format)
  {}

  void preprocess() {}
  void postprocess()
  {
    this->m_finish = true;
    if (this->m_clear)
    {
      for (auto& f : KmDir::get().get_files_to_merge(m_part_id, m_lz4, KM_FILE::KMER))
      {
        //std::remove(f.c_str());
        Eraser::get().erase(f);
      }
    }
    this->exec_callback();
  }

  void exec()
  {
    spdlog::debug("[exec] - KmerMergeTask - P={}", m_part_id);

    std::vector<std::string> paths = KmDir::get().get_files_to_merge(m_part_id,
                                                                     m_lz4,
                                                                     KM_FILE::KMER);
    std::string out_path = KmDir::get().get_matrix_path(m_part_id, m_mode, m_format,
                                                        COUNT_FORMAT::KMER, m_lz4);
    KmerMerger<span, MAX_C> merger(paths, m_ab_vec, m_kmer_size, m_rec_min, m_save_if);

#ifdef WITH_PLUGIN
    IMergePlugin* plugin = nullptr;

    if (PluginManager<IMergePlugin>::get().use_plugin())
    {
      plugin = PluginManager<IMergePlugin>::get().get_plugin();
      plugin->set_out_dir(KmDir::get().m_plugin_storage);
      plugin->set_kmer_size(m_kmer_size);
      plugin->set_partition(m_part_id);
      merger.set_plugin(plugin);
    }
#endif

    if (m_mode == MODE::COUNT)
    {
      if (m_format == FORMAT::TEXT)
        merger.write_as_text(out_path);
      else if (m_format == FORMAT::BIN)
        merger.write_as_bin(out_path, m_lz4);
    }
    else if (m_mode == MODE::PA)
    {
      if (m_format == FORMAT::TEXT)
        merger.write_as_pa_text(out_path);
      else if (m_format == FORMAT::BIN)
        merger.write_as_pa(out_path, m_lz4);
    }

#ifdef WITH_PLUGIN
    if (PluginManager<IMergePlugin>::get().use_plugin())
    {
      PluginManager<IMergePlugin>::get().destroy_plugin(plugin);
    }
    else
    {
      merger.get_infos()->serialize(KmDir::get().get_merge_info_path(m_part_id));
    }
#endif

    merger.get_infos()->serialize(KmDir::get().get_merge_info_path(m_part_id));

    spdlog::debug("[done] - KmerMergeTask - P={}", m_part_id);
  }

private:
  uint32_t m_part_id;
  std::vector<uint32_t>& m_ab_vec;
  uint32_t m_kmer_size;
  uint32_t m_rec_min;
  uint32_t m_save_if;
  bool m_lz4;
  MODE m_mode;
  FORMAT m_format;
};

template<size_t MAX_C>
class HashMergeTask : public ITask
{
public:
  HashMergeTask(uint32_t partition_id,
                std::vector<uint32_t>& ab_vec,
                uint32_t recurrence_min,
                uint32_t save_if,
                bool lz4,
                MODE mode,
                FORMAT format,
                HashWindow& win,
                bool clear,
                int32_t bw)
  : ITask(4, clear), m_part_id(partition_id), m_ab_vec(ab_vec), m_rec_min(recurrence_min),
    m_save_if(save_if), m_lz4(lz4), m_mode(mode), m_format(format), m_win(win), m_bw(bw) {}

  void preprocess() {}
  void postprocess()
  {
    if (this->m_clear)
    {
      for (auto& f : KmDir::get().get_files_to_merge(m_part_id, m_lz4, KM_FILE::HASH))
      {
        Eraser::get().erase(f);
      }
    }
    this->m_finish = true;
    this->exec_callback();
  }

  void exec()
  {
    spdlog::debug("[exec] - HashMergeTask - P={}", m_part_id);

    std::vector<std::string> paths = KmDir::get().get_files_to_merge(m_part_id,
                                                                     m_lz4,
                                                                     KM_FILE::HASH);
    std::string out_path = KmDir::get().get_matrix_path(m_part_id, m_mode, m_format,
                                                        COUNT_FORMAT::HASH, false);

    HashMerger<MAX_C, 32768, HashReader<MAX_C, 32768>> merger(paths, m_ab_vec, m_rec_min, m_save_if);

#ifdef WITH_PLUGIN
    IMergePlugin* plugin = nullptr;

    if (PluginManager<IMergePlugin>::get().use_plugin())
    {
      plugin = PluginManager<IMergePlugin>::get().get_plugin();
      plugin->set_out_dir(KmDir::get().m_plugin_storage);
      plugin->set_kmer_size(0);
      plugin->set_partition(m_part_id);
      merger.set_plugin(plugin);
    }
#endif

    if (m_mode == MODE::COUNT)
    {
      if (m_format == FORMAT::TEXT)
        merger.write_as_text(out_path);
      else if (m_format == FORMAT::BIN)
        merger.write_as_bin(out_path, m_lz4);
    }
    else if (m_mode == MODE::PA)
    {
      if (m_format == FORMAT::TEXT)
        merger.write_as_pa_text(out_path);
      else if (m_format == FORMAT::BIN)
        merger.write_as_pa(out_path, m_lz4);
    }
    else if (m_mode == MODE::BF)
    {
        merger.write_as_bf(out_path, m_win.get_lower(m_part_id),
                           m_win.get_upper(m_part_id), false);
    }
    else if (m_mode == MODE::BFT)
    {
        merger.write_as_bft(out_path, m_win.get_lower(m_part_id),
                            m_win.get_upper(m_part_id), false);
    }
    else if (m_mode == MODE::BFC)
    {
        merger.write_as_bfc(out_path, m_win.get_lower(m_part_id),
                            m_win.get_upper(m_part_id), m_bw, false);
    }

#ifdef WITH_PLUGIN
    if (PluginManager<IMergePlugin>::get().use_plugin())
    {
      PluginManager<IMergePlugin>::get().destroy_plugin(plugin);
    }
    else
    {
      merger.get_infos()->serialize(KmDir::get().get_merge_info_path(m_part_id));
    }
#endif
    merger.get_infos()->serialize(KmDir::get().get_merge_info_path(m_part_id));

    if (m_mode == MODE::BF || m_mode == MODE::BFT)
    {
      std::string fpr_path = fmt::format("{}/{}", KmDir::get().m_fpr_storage, fmt::format("partition_{}.txt", m_part_id));
      std::ofstream fp(fpr_path, std::ios::out); check_fstream_good(fpr_path, fp);

      size_t m = m_win.get_window_size_bits();
      for (auto& n : merger.get_infos()->get_unique_w_rescue())
      {
        double fpr = bloom_fp(m, n);
        fp << std::fixed << fpr << "\n";
      }
    }

    spdlog::debug("[done] - HashMergeTask - P={}", m_part_id);
  }

private:
  uint32_t m_part_id;
  std::vector<uint32_t>& m_ab_vec;
  uint32_t m_rec_min;
  uint32_t m_save_if;
  bool m_lz4;
  MODE m_mode;
  FORMAT m_format;
  HashWindow& m_win;
  uint32_t m_bw;
};


class FormatVectorTask : public ITask
{
public:
  FormatVectorTask(std::string id, OUT_FORMAT bf_type, uint64_t bloom,
                   uint32_t nb_parts, bool lz4, uint32_t kmer_size, bool clear = false)
    : ITask(5, clear),
      m_id(id), m_bf_type(bf_type), m_nb_parts(nb_parts), m_lz4(lz4), m_bloom(bloom),
      m_kmer_size(kmer_size)
  {}

  void preprocess() {}
  void postprocess()
  {
    if (this->m_clear)
    {
      for (size_t p=0; p<m_nb_parts; p++)
      {
        std::string s = KmDir::get().get_count_part_path(m_id, p, this->m_lz4, KM_FILE::VECTOR);
        Eraser::get().erase(s);
      }
    }
    this->m_finish = true;
    this->exec_callback();
  }

  void exec()
  {
    spdlog::debug("[exec] - FormatVectorTask - S={}", m_id);
    BloomBuilderFromVec(KmDir::get().m_fof.get_i(m_id), m_bf_type, m_bloom, m_nb_parts, m_kmer_size, m_lz4).build();
    spdlog::debug("[done] - FormatVectorTask - S={}", m_id);
  }

private:
  std::string m_id;
  OUT_FORMAT m_bf_type;
  uint32_t m_nb_parts;
  bool m_lz4;
  uint64_t m_bloom;
  uint32_t m_kmer_size;
};

class FormatTask : public ITask
{
public:
  FormatTask(std::vector<int>& files, std::vector<std::mutex>& file_mutex,
              OUT_FORMAT bf_type, uint64_t bloom, uint32_t file_id, uint32_t nb_parts,
              uint32_t kmer_size, bool clear = false)
    : ITask(5, clear), m_fds(files), m_mutex(file_mutex),
      m_bf_type(bf_type), m_file_id(file_id), m_nb_parts(nb_parts), m_bloom(bloom),
      m_kmer_size(kmer_size)
  {}

  void preprocess() {}
  void postprocess()
  {
    this->m_finish = true;
    this->exec_callback();
  }

  void exec()
  {
    spdlog::debug("[exec] - FormatTask - S={}", KmDir::get().m_fof.get_id(m_file_id));
    BloomBuilderFromHash(m_fds, m_mutex, m_bf_type, m_bloom, m_file_id, m_nb_parts, m_kmer_size).build();
    spdlog::debug("[done] - FormatTask - S={}", KmDir::get().m_fof.get_id(m_file_id));
  }

private:
  OUT_FORMAT m_bf_type;
  uint32_t m_nb_parts;
  uint32_t m_file_id;
  uint64_t m_bloom;
  uint32_t m_kmer_size;
  std::vector<int> m_fds;
  std::vector<std::mutex>& m_mutex;
};

};
