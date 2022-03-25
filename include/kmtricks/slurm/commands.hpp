#pragma once

#include <string>
#include <filesystem>
#include <kmtricks/io/fof.hpp>
#include <kmtricks/cmd/all.hpp>

#include <fmt/format.h>

namespace km {

  namespace slurm {

    namespace fs = std::filesystem;

    const std::string repart_cmd =
      "srun {} repart --file {} --run-dir {} --kmer-size {} --minimizer-size {} "
      "--nb-partitions {} --bloom-size {} --verbose {}";

    const std::string superk_cmd =
      "srun {} superk --run-dir {} --id {} --verbose {}";

    const std::string count_cmd =
      "srun {} count --run-dir {} --id {} --hard-min {} --mode {} --threads {} --verbose {}";

    const std::string merge_cmd =
      "srun {} merge --run-dir {} --soft-min {} --recurrence-min {} --share-min {} "
      " --mode {} --threads {} --verbose {}";

    const std::string format_cmd =
      "srun {} format --run-dir {} --out-format {} --threads {} --verbose {}";

    const std::string submit_cmd =
      "{}=$(sbatch {})";

    const std::string submit_cmd_d =
      "{}=$(sbatch --dependency=afterok:${{{}##* }} {})";

    const std::string submit_cmd_d_m =
      "echo $(sbatch --array=$(ls {} | echo \"$(wc -l)-1\" | bc)%{} {})";

    const std::string submit_cmd_d_f =
      "sbatch --dependency=afterok:$(cat {} | cut -d' ' -f4) {}";

    const std::string r_script = "{}/kmtricks_repart.slurm";
    const std::string sk_script = "{}/kmtricks_sk.slurm";
    const std::string skc_script = "{}/kmtricks_skc.slurm";
    const std::string m_script = "{}/kmtricks_merge.slurm";
    const std::string f_script = "{}/kmtricks_format.slurm";
    const std::string submit_script = "{}/submit.sh";
    const std::string submit_script_m = "{}/submit_merge.slurm";
    const std::string submit_script_f = "{}/submit_format.slurm";

    const std::string slurm_array = "#SBATCH --array=0-{}%{}\n";

    const std::string shebang = "#!/bin/bash";

    const std::string slurm_template =
      "#SBATCH --job-name={}\n"
      "#SBATCH --cpus-per-task={}\n"
      "#SBATCH --mem-per-cpu={}\n"
      "#SBATCH --ntasks={}\n"
      "#SBATCH --error={}\n"
      "#SBATCH --output={}\n";

    void add_slurm_options(std::ostream& s, const std::vector<std::pair<std::string, std::string>>& options)
    {
      for (auto& [k, v] : options)
        s << fmt::format("#SBATCH {}={}\n", k, v);
      s << std::endl;
    }

    std::string exec_path()
    {
      char path[512] = {0};
      ::readlink("/proc/self/exe", path, 512);
      return std::string(path);
    }

    std::string abspath(const std::string& p)
    {
      return fs::absolute(fs::path(p));
    }

    const std::string kmtricks_bin = exec_path();

    std::string id_from_fof(const std::string& fof_path)
    {
      return fmt::format("$(awk -v sample={} 'NR==sample+1' {} | cut -d':' -f1)", "${SLURM_ARRAY_TASK_ID}", fof_path);
    }

    void repart(const std::string& path,
                       std::vector<std::string>& cmds,
                       all_options_t opt)
    {
      std::ofstream out(path, std::ios::out);

      out << shebang << "\n\n";

      out << fmt::format(slurm_template, "kmr",
                                         1,
                                         "5G",
                                         1,
                                         fmt::format("{}/%x_%j.err", opt->slurm_dir),
                                         fmt::format("{}/%x_%j.out", opt->slurm_dir));

      add_slurm_options(out, opt->slurm_options);

      out << std::endl;

      out << fmt::format(repart_cmd, kmtricks_bin,
                                     opt->fof,
                                     opt->dir,
                                     opt->kmer_size,
                                     opt->minim_size,
                                     opt->nb_parts,
                                     opt->bloom_size,
                                     opt->verbosity);

      cmds.push_back(fmt::format(submit_cmd, "KMR_ID", fmt::format(r_script, opt->slurm_dir)));
    }

    void superk(const std::string& path,
                       std::vector<std::string>& cmds,
                       all_options_t opt)
    {
      std::ofstream out(path, std::ios::out);

      out << shebang << "\n\n";
      out << fmt::format(slurm_template, "kmsk",
                                         4,
                                         "5G",
                                         1,
                                         fmt::format("{}/%x_%A_%a.err", opt->slurm_dir),
                                         fmt::format("{}/%x_%A_%a.out", opt->slurm_dir));

      out << fmt::format(slurm_array, Fof(opt->fof).size() - 1, opt->slurm_max_array);
      add_slurm_options(out, opt->slurm_options);
      out << fmt::format(superk_cmd, kmtricks_bin,
                                     opt->dir,
                                     id_from_fof(abspath(opt->fof)),
                                     opt->verbosity);

      out << (opt->lz4 ? " --cpr" : "");

      out << std::endl;

      cmds.push_back(fmt::format(submit_cmd_d, "KMSK_ID", "KMR_ID", fmt::format(sk_script, opt->slurm_dir)));
    }

    void superk_count(const std::string& path,
                             std::vector<std::string>& cmds,
                             all_options_t opt)
    {
      std::ofstream out(path, std::ios::out);
      out << shebang << "\n\n";
      out << fmt::format(slurm_template, "kmskc",
                                         opt->nb_threads,
                                         opt->slurm_mem,
                                         1,
                                         fmt::format("{}/%x_%A_%a.err", opt->slurm_dir),
                                         fmt::format("{}/%x_%A_%a.out", opt->slurm_dir));

      out << fmt::format(slurm_array, Fof(opt->fof).size() - 1, opt->slurm_max_array);
      add_slurm_options(out, opt->slurm_options);

      out << fmt::format(superk_cmd, kmtricks_bin,
                                     opt->dir,
                                     id_from_fof(abspath(opt->fof)),
                                     opt->verbosity);
      out << (opt->lz4 ? " --cpr" : "");

      out << std::endl;

      std::string mode = opt->skip_merge ? "vector" : cformat_to_str(opt->count_format);

      out << fmt::format(count_cmd, kmtricks_bin,
                                    opt->dir,
                                    id_from_fof(abspath(opt->fof)),
                                    opt->c_ab_min,
                                    mode,
                                    opt->nb_threads,
                                    opt->verbosity);

      out << (opt->hist ? " --hist" : "");
      out << (opt->lz4 ? " --cpr" : "");
      out << (opt->keep_tmp ? "" : " --clear");

      out << std::endl;

      cmds.push_back(fmt::format(submit_cmd_d, "KMSKC_ID", "KMR_ID", fmt::format(skc_script, opt->slurm_dir)));
    }

    void merge(const std::string& path,
                      std::vector<std::string>& cmds,
                      all_options_t opt)
    {
      std::ofstream out(path, std::ios::out);
      out << shebang << "\n\n";
      out << fmt::format(slurm_template, "kmm",
                                         opt->nb_threads,
                                         "500M",
                                         1,
                                         fmt::format("{}/%x_%A_%a.err", opt->slurm_dir),
                                         fmt::format("{}/%x_%A_%a.out", opt->slurm_dir));

      add_slurm_options(out, opt->slurm_options);

      std::string mode = fmt::format("{}:{}:{}", cformat_to_str(opt->count_format), mode_to_str(opt->mode), format_to_str2(opt->format));

      std::string soft_min;
      if (!opt->m_ab_min_path.empty())
        soft_min = opt->m_ab_min_path;
      else if (opt->m_ab_float)
        soft_min = std::to_string(opt->m_ab_min_f);
      else
        soft_min = std::to_string(opt->m_ab_min);

      out << fmt::format(merge_cmd, kmtricks_bin,
                                    opt->dir,
                                    opt->m_ab_min,
                                    opt->r_min,
                                    soft_min,
                                    mode,
                                    opt->nb_threads,
                                    opt->verbosity);
      out << (opt->lz4 ? " --cpr" : "");
      out << (opt->keep_tmp ? "" : " --clear");

      out << std::endl;

      std::ofstream out2(fmt::format(submit_script_m, opt->slurm_dir), std::ios::out);

      out2 << shebang << "\n\n";

      out2 << "#SBATCH --job-name=kmm_submit\n";
      out2 << fmt::format("#SBATCH --output={}/{}\n\n", opt->slurm_dir, "MERGE_PID");

      out2 << fmt::format(submit_cmd_d_m,
                          fmt::format("{}/minimizers", opt->dir),
                          opt->slurm_max_array,
                          fmt::format(m_script, opt->slurm_dir));
      out2 << std::endl;

      cmds.push_back(fmt::format(submit_cmd_d,
                                 "KMM_ID",
                                 "KMSKC_ID",
                                 fmt::format(submit_script_m, opt->slurm_dir)));
    }

    void format_bf(const std::string& path,
                          std::vector<std::string>& cmds,
                          all_options_t opt)
    {
      std::ofstream out(path, std::ios::out);
      out << shebang << "\n\n";
      out << fmt::format(slurm_template, "kmf",
                                         opt->nb_threads,
                                         "500M",
                                         1,
                                         fmt::format("{}/%x_%j.err", opt->slurm_dir),
                                         fmt::format("{}/%x_%j.out", opt->slurm_dir)) << "\n\n";

      add_slurm_options(out, opt->slurm_options);
      out << fmt::format(format_cmd, kmtricks_bin,
                                     opt->dir,
                                     format_to_str(opt->out_format),
                                     opt->nb_threads,
                                     opt->verbosity);

      out << (opt->skip_merge ? " --from-vec" : " --from-hash");
      out << (opt->lz4 ? " --cpr-in" : "");
      out << (opt->keep_tmp ? "" : " --clear");

      out << std::endl;

      std::ofstream out2(fmt::format(submit_script_f, opt->slurm_dir), std::ios::out);
      out2 << shebang << "\n\n";
      out2 << "#SBATCH --job-name=kmf_submit\n";
      out2 << fmt::format(submit_cmd_d_f, fmt::format("{}/MERGE_PID", opt->slurm_dir), fmt::format(f_script, opt->slurm_dir));

      cmds.push_back(fmt::format(submit_cmd_d, "KMF", "KMM_ID", fmt::format(submit_script_f, opt->slurm_dir)));
    }

    void submit(const std::string& path, const std::vector<std::string>& cmds)
    {
      std::ofstream out(path, std::ios::out);
      out << shebang << "\n\n";

      for (auto& c : cmds)
      {
        out << c << "\n";
      }
    }

  } // end of namespace slurm
} // end of namespace km
