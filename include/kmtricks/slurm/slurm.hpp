#pragma once

#include <string>
#include <vector>

#include <kmtricks/slurm/commands.hpp>
#include <kmtricks/slurm/options.hpp>

namespace km {
  namespace slurm {

    void slurm(all_options_t opt)
    {
      spdlog::warn("kmtricks slurm support is currently experimental.");

      opt->dir = abspath(opt->dir);
      opt->fof = abspath(opt->fof);

      fs::create_directory(fs::path(opt->slurm_dir));

      opt->slurm_dir = abspath(opt->slurm_dir);

      std::vector<std::string> cmds;

      valid_slurm_options(opt->slurm_options);

      repart(fmt::format(r_script, opt->slurm_dir), cmds, opt);

      if (opt->until == COMMAND::SUPERK)
      {
        superk(fmt::format(sk_script, opt->slurm_dir), cmds, opt);
        return;
      }

      if (opt->until == COMMAND::COUNT || opt->until == COMMAND::ALL)
        superk_count(fmt::format(skc_script, opt->slurm_dir), cmds, opt);

      if (opt->until == COMMAND::MERGE || opt->until == COMMAND::ALL)
        merge(fmt::format(m_script, opt->slurm_dir), cmds, opt);

      if ((opt->until == COMMAND::FORMAT || opt->until == COMMAND::ALL) && opt->count_format == COUNT_FORMAT::HASH)
        format_bf(fmt::format(f_script, opt->slurm_dir), cmds, opt);

      submit(fmt::format(submit_script, opt->slurm_dir), cmds);

      if (opt->slurm_submit)
      {
        std::system(fmt::format("bash {}/submit.sh &", opt->slurm_dir).c_str());
        spdlog::info("Submitted. (logs at {})", opt->slurm_dir);
      }
      else
      {
        spdlog::info("Done. To submit: `bash {}/submit.sh`", opt->slurm_dir);
      }
    }

  } // end of namespace slurm
} // end of namespace km
