#pragma once

#include <set>

namespace km {
  namespace slurm {

    const std::set<std::string> slurm_options_set {
      "array", "a",
      "account", "A",
      "acctg_freq",
      "extra_node_info", "B",
      "batch",
      "bb",
      "bbf",
      "begin", "b",
      "cluster_constraint",
      "comment",
      "constraint", "C",
      "contiguous",
      "cores_per_socket",
      "cpu_freq",
      "cpus_per_gpu",
      "cpus_per_task", "c",
      "deadline"
      "delay_boot"
      "dependency", "d",
      "chdir", "D",
      "error", "e",
      "exclusive",
      "export",
      "export_file",
      "nodefile", "F",
      "get_user_env",
      "gid",
      "gpus", "G"
      "gpu_bind",
      "gpu_freq",
      "gpus_per_node",
      "gpus_per_socket",
      "gpus_per_task",
      "gres",
      "gres_flags",
      "hold", "H",
      "hint",
      "ignore_pbs",
      "input", "i",
      "job_name", "J",
      "no_kill", "k",
      "kill_on_invalid_dep",
      "licenses", "L",
      "clusters", "M",
      "distribution", "m",
      "mail_type",
      "mail_user",
      "mcs_label",
      "mem",
      "mem_per_cpu",
      "mem_per_gpu",
      "mem_bind",
      "mincpus",
      "nodes", "N",
      "ntasks", "n",
      "network",
      "nice",
      "no_requeue",
      "ntasks_per_core",
      "ntasks_per_node",
      "ntasks_per_socket",
      "overcommit", "O",
      "output", "o",
      "open_mode",
      "parsable",
      "partition", "p",
      "power",
      "priority",
      "profile",
      "propagate",
      "qos", "q",
      "quiet", "Q",
      "reboot",
      "requeue",
      "reservation",
      "oversubscribe", "s",
      "core_spec", "S",
      "signal",
      "sockets_per_node",
      "spread_job",
      "switches",
      "time", "t",
      "test_only",
      "thread_spec",
      "threads_per_core",
      "time_min",
      "tmp",
      "usage",
      "uid",
      "use_min_nodes",
      "version", "V",
      "verbose", "v",
      "nodelist", "w",
      "wait", "W",
      "wait_all_nodes",
      "wckey",
      "wrap",
      "exclude", "x"
    };

    const std::set<std::string> slurm_reserved_options {
      "cpus-per-task",
      "mem-per-cpu",
      "job-name",
      "ntasks",
      "error",
      "output"
    };

    bool is_valid_slurm_opt(const std::string& option)
    {
      return static_cast<bool>(slurm_options_set.count(option));
    }

    bool is_free_slurm_opt(const std::string& option)
    {
      return !static_cast<bool>(slurm_reserved_options.count(option));
    }

    void add_prefix(std::string& o)
    {
      if (o.size() > 1)
        o = fmt::format("--{}", o);
      else
        o = fmt::format("-{}", o);
    }

    km_EXCEPTION(SlurmError);

    void valid_slurm_options(std::vector<std::pair<std::string, std::string>>& options)
    {
      for (auto& [o, v] : options)
      {
        if (!is_valid_slurm_opt(o))
          throw SlurmError(fmt::format("Invalid slurm option: '{}={}'.", o, v));
        else if (!is_free_slurm_opt(o))
          throw SlurmError(fmt::format("'{}' is set internally by kmtricks.", o));
        else
          add_prefix(o);
      }
    }


  } // end of namespace slurm
} // end of namespace km
