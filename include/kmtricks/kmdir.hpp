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
#include <filesystem>

#include <kmtricks/cmd/cmd_common.hpp>
#include <kmtricks/cmd/infos.hpp>
#include <kmtricks/io/io_common.hpp>
#include <kmtricks/io/fof.hpp>
#include <fmt/format.h>

namespace fs = std::filesystem;

namespace km {

class KmDir
{
public:
  static KmDir& get()
  {
    static KmDir kmdir;
    return kmdir;
  }

  std::string get_superk_path(const std::string& sample_id)
  {
    return fmt::format("{}/{}", m_superk_storage, sample_id);
  }

  std::vector<std::string> get_files_to_merge(uint32_t part_id, bool compressed, KM_FILE km_file)
  {
    std::string ext;

    if (KM_FILE::HASH == km_file)
    {
      ext = "hash";
      if (compressed) ext += ".p4";
    }
    else if (KM_FILE::KMER == km_file)
    {
      ext = "kmer";
      if (compressed) ext += ".lz4";
    }

    std::vector<std::string> paths;

    Fof fof(m_fof_path);
    for (auto s : fof)
    {
      paths.push_back(fmt::format(m_part_template, m_counts_storage, part_id, std::get<0>(s), ext));
      if (!fs::exists(paths.front()))
        throw FileNotFoundError(fmt::format("{} is missing.", paths.front()));
    }
    return paths;
  }

  std::string get_count_part_path(std::string id,
                                  uint32_t part_id,
                                  bool compressed,
                                  KM_FILE km_file)
  {
    std::string ext;
    if (KM_FILE::HASH == km_file)
      ext = "hash";
    else if (KM_FILE::KMER == km_file)
      ext = "kmer";
    else if (KM_FILE::VECTOR == km_file)
      ext = "vector";
    else if (KM_FILE::KFF == km_file)
      ext = "kff";

    if (compressed)
    {
      if (KM_FILE::KMER == km_file || KM_FILE::VECTOR == km_file)
        ext += ".lz4";
      else if (KM_FILE::HASH == km_file)
        ext += ".p4";
    }

    return fmt::format(m_part_template, m_counts_storage, part_id, id, ext);
  }

  std::vector<std::string> get_count_part_paths(std::string id, uint32_t nb_parts,
                                                bool compressed, KM_FILE km_file)
  {
    std::vector<std::string> paths;
    for (size_t i=0; i<nb_parts; i++)
    {
      std::string p = get_count_part_path(id, i, compressed, km_file);
      if (fs::exists(p))
        paths.push_back(p);
    }
    return paths;
  }

  std::string get_matrix_path(uint32_t part_id, MODE mode, FORMAT format,
                              COUNT_FORMAT cformat, bool compressed)
  {
    std::string ext;
    if (MODE::COUNT == mode)
      ext = cformat == COUNT_FORMAT::KMER ? "count" : "count_hash";
    else if (MODE::PA == mode)
      ext = cformat == COUNT_FORMAT::KMER ? "pa" : "pa_hash";
    else if (MODE::BF == mode)
      ext = "cmbf";
    else if (MODE::BFT == mode)
      ext = "rmbf";
    else if (MODE::BFC == mode)
      ext = "cmbf";

    if (FORMAT::TEXT == format)
      ext += ".txt";

    if (compressed && (mode != MODE::BFT) && format != FORMAT::TEXT)
      ext += ".lz4";

    return fmt::format(m_matrix_template, m_matrix_storage, part_id, ext);
  }

  std::vector<std::string> get_matrix_paths(uint32_t nb_parts, MODE mode,
                                           FORMAT format, COUNT_FORMAT cformat, bool compressed)
  {
    std::vector<std::string> paths;
    for (size_t i=0; i<nb_parts; i++)
    {
      std::string p = get_matrix_path(i, mode, format, cformat, compressed);
      if (fs::exists(p))
        paths.push_back(p);
    }
    return paths;
  }

  std::string get_filter_path(const std::string& id, OUT_FORMAT out)
  {
    return fmt::format(m_filter_template, m_filter_storage, id,
                       out == OUT_FORMAT::HOWDE ? "bf" : "sdsl");
  }


  std::string get_hist_path(const std::string& id)
  {
    return fmt::format(m_hist_template, m_hist_storage, id);
  }

  std::string get_merge_info_path(uint32_t part_id)
  {
    return fmt::format(m_stat_merge_template, m_stat_storage, part_id);
  }

  std::string get_bf_list_path()
  {
    return fs::absolute(fs::path(fmt::format("{}/bf_list", m_index_storage))).string();
  }

  std::string get_index_path()
  {
    return fs::absolute(fs::path(fmt::format("{}/index", m_index_storage))).string();
  }

  std::string get_pinfos_path(const std::string& id)
  {
    return fmt::format("{}/{}.pinfo", m_part_info_storage, id);
  }

  std::string get_merge_th_path()
  {
    return fmt::format("{}/merge_amin.txt", m_root);
  }

  std::vector<std::string> get_minim_paths(uint32_t nb_parts)
  {
    fs::create_directory(m_minimizer_storage);
    std::vector<std::string> paths;
    for (size_t i=0; i<nb_parts; i++)
      paths.push_back(fmt::format("{}/minimizers.{}", m_minimizer_storage, i));
    return paths;
  }

  void init(const std::string& root, const std::string& fof, bool first = false)
  {
    m_root = fs::absolute(fs::path(root)).string();
    m_fof_path = fmt::format("{}/kmtricks.fof", m_root);
    m_config_storage = fmt::format("{}/config", m_root);
    m_repart_storage = fmt::format("{}/repartition", m_root);
    m_superk_storage = fmt::format("{}/superkmers", m_root);
    m_counts_storage = fmt::format("{}/counts", m_root);
    m_matrix_storage = fmt::format("{}/matrices", m_root);
    m_filter_storage = fmt::format("{}/filters", m_root);
    m_hist_storage = fmt::format("{}/histograms", m_root);
    m_stat_storage = fmt::format("{}/merge_infos", m_root);
    m_index_storage = fmt::format("{}/howde_index", m_root);
    m_part_info_storage = fmt::format("{}/partition_infos", m_root);
    m_hash_win = fmt::format("{}/hash.info", m_root);
    m_run_infos = fmt::format("{}/run_infos.txt", m_root);
    m_options = fmt::format("{}/options.txt", m_root);
    m_minimizer_storage = fmt::format("{}/minimizers", m_root);
    m_fpr_storage = fmt::format("{}/fpr", m_root);
    m_plugin_storage = fmt::format("{}/plugin_output", m_root);

    if (first)
    {
      m_fof = Fof(fof);
      fs::create_directory(m_root);
      m_fof.copy(m_fof_path);
      fs::create_directory(m_superk_storage);
      fs::create_directory(m_counts_storage);
      fs::create_directory(m_matrix_storage);
      fs::create_directory(m_filter_storage);
      fs::create_directory(m_hist_storage);
      fs::create_directory(m_stat_storage);
      fs::create_directory(m_index_storage);
      fs::create_directory(m_part_info_storage);
      fs::create_directory(m_fpr_storage);
#ifdef WITH_PLUGIN
      fs::create_directory(m_plugin_storage);
#endif
      std::string info_path = fmt::format("{}/build_infos.txt", m_root);
      std::ofstream bout(info_path, std::ios::out); check_fstream_good(info_path, bout);
      main_infos(bout);
    }
    else
    {
      m_fof = Fof(m_fof_path);
    }
  }

  void init_part(uint32_t nb_parts)
  {
    for (uint32_t i=0; i<nb_parts; i++)
    {
      fs::create_directory(fmt::format("{}/partition_{}", m_counts_storage, i));
    }
  }

  void init_one_part(uint32_t part)
  {
    fs::create_directory(fmt::format("{}/partition_{}", m_counts_storage, part));
  }

private:
  KmDir() {}

public:
  std::string m_root;
  std::string m_fof_path;
  std::string m_config_storage;
  std::string m_repart_storage;
  std::string m_superk_storage;
  std::string m_counts_storage;
  std::string m_matrix_storage;
  std::string m_filter_storage;
  std::string m_hist_storage;
  std::string m_stat_storage;
  std::string m_index_storage;
  std::string m_hash_win;
  std::string m_part_info_storage;
  std::string m_minimizer_storage;
  std::string m_run_infos;
  std::string m_options;
  std::string m_fpr_storage;
  std::string m_plugin_storage;

  std::string m_filter_template {"{}/{}.{}"};
  std::string m_matrix_template {"{}/matrix_{}.{}"};
  std::string m_part_template {"{}/partition_{}/{}.{}"};
  std::string m_hist_template {"{}/{}.hist"};
  std::string m_stat_merge_template {"{}/partition{}.merge_info"};

  Fof m_fof;
};

};
