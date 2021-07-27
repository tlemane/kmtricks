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

#include <kmtricks/cli.hpp>
#include <kmtricks/utils.hpp>
#include <filesystem>

namespace fs = std::filesystem;
namespace km {

kmtricksCli::kmtricksCli(
    const std::string& name,
    const std::string& desc,
    const std::string& version,
    const std::string& authors)
{
  cli = std::make_shared<bc::Parser<1>>(bc::Parser<1>(name, desc, version, authors));
  all_opt = std::make_shared<struct all_options>(all_options{});
  repart_opt = std::make_shared<struct repart_options>(repart_options{});
  superk_opt = std::make_shared<struct superk_options>(superk_options{});
  count_opt = std::make_shared<struct count_options>(count_options{});
  merge_opt = std::make_shared<struct merge_options>(merge_options{});
  format_opt = std::make_shared<struct format_options>(format_options{});
  dump_opt = std::make_shared<struct dump_options>(dump_options{});
  agg_opt = std::make_shared<struct agg_options>(agg_options{});
  index_opt = std::make_shared<struct index_options>(index_options{});
  query_opt = std::make_shared<struct query_options>(query_options{});
  all_cli(cli, all_opt);
#ifdef WITH_KM_MODULES
  repart_cli(cli, repart_opt);
  superk_cli(cli, superk_opt);
  count_cli(cli, count_opt);
  merge_cli(cli, merge_opt);
  format_cli(cli, format_opt);
#endif
  dump_cli(cli, dump_opt);
  agg_cli(cli, agg_opt);
#ifdef WITH_HOWDE
  index_cli(cli, index_opt);
  query_cli(cli, query_opt);
#endif
  info_cli(cli);
}

std::tuple<COMMAND, km_options_t> kmtricksCli::parse(int argc, char* argv[])
{
  if (argc > 1 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h"))
  {
    cli->show_help();
    exit(EXIT_FAILURE);
  }
  else if (argc > 1 && (std::string(argv[1]) == "--version" || std::string(argv[1]) == "-v"))
  {
    std::cerr << fmt::format("kmtricks {}", PROJECT_VER) << std::endl;
    exit(EXIT_FAILURE);
  }

  try
  {
    (*cli).parse(argc, argv);
  }
  catch (const bc::ex::BCliError& e)
  {
    bc::utils::exit_bcli(e);
    exit(EXIT_FAILURE);
  }

  if (cli->is("pipeline"))
    return std::make_tuple(COMMAND::ALL, all_opt);
  else if (cli->is("repart"))
    return std::make_tuple(COMMAND::REPART, repart_opt);
  else if (cli->is("superk"))
    return std::make_tuple(COMMAND::SUPERK, superk_opt);
  else if (cli->is("count"))
    return std::make_tuple(COMMAND::COUNT, count_opt);
  else if (cli->is("merge"))
    return std::make_tuple(COMMAND::MERGE, merge_opt);
  else if (cli->is("format"))
    return std::make_tuple(COMMAND::FORMAT, format_opt);
  else if (cli->is("dump"))
    return std::make_tuple(COMMAND::DUMP, dump_opt);
  else if (cli->is("aggregate"))
    return std::make_tuple(COMMAND::AGGREGATE, agg_opt);
  else if (cli->is("index"))
    return std::make_tuple(COMMAND::INDEX, index_opt);
  else if (cli->is("query"))
    return std::make_tuple(COMMAND::QUERY, query_opt);
  else
    return std::make_tuple(COMMAND::INFOS, std::make_shared<struct km_options>(km_options{}));
}

auto dir_already_exists = [](const std::string& p, const std::string& v) {
  bool exists = !fs::is_directory(v);
  return std::make_tuple(exists, bc::utils::format_error(p, v, "Directory already exists!"));
};

km_options_t all_cli(std::shared_ptr<bc::Parser<1>> cli, all_options_t options)
{
  bc::cmd_t all_cmd = cli->add_command("pipeline", "kmtricks pipeline (run all the steps, repart -> superk -> count -> merge -> format)");

  all_cmd->add_param("--file", "kmtricks input file, see README.md.")
    ->meta("FILE")
    ->checker(bc::check::is_file)
    ->setter(options->fof);

  all_cmd->add_param("--run-dir", "kmtricks runtime directory.")
    ->meta("DIR")
    ->checker(dir_already_exists)
    ->setter(options->dir);

  all_cmd->add_param("--kmer-size", fmt::format("size of a k-mer. [8, {}].", KL[KMER_N-1]-1))
    ->meta("INT")
    ->def("31")
    ->checker(bc::check::f::range(8, KL[KMER_N-1]-1))
    ->setter(options->kmer_size);

  all_cmd->add_param("--count-abundance-min", "min abundance to keep a k-mer.")
    ->meta("INT")
    ->def("2")
    ->checker(bc::check::is_number)
    ->setter(options->c_ab_min);

  auto mode_setter = [options](const std::string& v) {
    auto s = bc::utils::split(v, ':');
    options->count_format = str_to_cformat(s[0]);
    options->mode = str_to_mode(s[1]);
    options->format = str_to_format2(s[2]);
  };

  auto mode_checker = [](const std::string& p, const std::string& v) -> bc::check::checker_ret_t {
    std::string available = "kmer:pa:text|"
                            "kmer:pa:bin|"
                            "kmer:count:text|"
                            "kmer:count:bin|"
                            "hash:count:text|"
                            "hash:count:bin|"
                            "hash:pa:text|"
                            "hash:pa:bin|"
                            "hash:bf:bin|"
                            "hash:bft:bin";
    auto s = bc::utils::split(v, ':');
    std::string mode;
    std::string format;
    std::string out;

    if (s.size() != 3)
      goto fail;

    mode = s[0];
    format = s[1];
    out = s[2];

    if (format != "count" && format != "pa" && format != "bf" && format != "bft")
      goto fail;

    if (mode == "kmer")
    {
      if (format == "bf" || format == "bft")
        goto fail;
    }
    else if (mode != "hash")
      goto fail;

    if ((format == "bf" || format == "bft") && out == "text")
      goto fail;

    goto success;

    fail:
      return std::make_tuple(false, bc::utils::format_error(
        p, v, fmt::format("Available -> {}", available)));

    success:
      return std::make_tuple(true, "");
  };

  all_cmd->add_param("--mode", "matrix mode <mode:format:out>, see README")
    ->meta("MODE:FORMAT:OUT")
    ->def("kmer:count:bin")
    ->checker(mode_checker)
    ->setter_c(mode_setter);

  all_cmd->add_param("--hist", "compute k-mer histograms.")
    ->as_flag()
    ->setter(options->hist);

  all_cmd->add_param("--kff-output", "output counted k-mers in kff format (only with --until count).")
    ->as_flag()
    ->setter(options->kff);

  all_cmd->add_param("--kff-sk-output", "output counted k-mers in kff-sk format (only with --until count).")
    ->as_flag()
    ->setter(options->kff_sk);

  all_cmd->add_param("--keep-tmp", "keep tmp files.")
    ->as_flag()
    ->setter(options->keep_tmp);

  all_cmd->add_group("merge options", "");

  auto a_min_setter = [options](const std::string& v) {
    auto [f, _] = bc::check::is_file("", v);
    if (f)
    {
      options->m_ab_min_path = v;
      return;
    }
    if (v.find('.') != std::string::npos)
    {
      bc::check::throw_if_false(bc::check::f::range(0.0, 1.0)("--abundance-min<float>", v));
      options->m_ab_min_f;
      options->m_ab_float = true;
    }
    else
      options->m_ab_min = bc::utils::lexical_cast<uint32_t>(v);
  };

  all_cmd->add_param("--merge-abundance-min", "during merge, min abundance to keep a k-mer.")
    ->meta("INT")
    ->def("1")
    ->setter_c(a_min_setter);

  all_cmd->add_param("--recurrence-min", "min recurrence to keep a k-mer.")
    ->meta("INT")
    ->def("1")
    ->checker(bc::check::is_number)
    ->setter(options->r_min);

  all_cmd->add_param("--save-if", "save a non-solid k-mer if it is solid in N other samples.")
    ->meta("INT")
    ->def("0")
    ->checker(bc::check::is_number)
    ->setter(options->save_if);

  all_cmd->add_group("pipeline control", "");

  auto until_setter = [options](const std::string& v) {
    options->until = str_to_cmd(v);
  };

  all_cmd->add_param("--until", "run until [all|repart|superk|count|merge|format]")
    ->meta("STR")
    ->def("all")
    ->checker(bc::check::f::in("all|repart|superk|count|merge|format"))
    ->setter_c(until_setter);

  all_cmd->add_param("--skip-merge", "skip merge step, only with --mode hash:bft:bin.")
    ->as_flag()
    ->setter(options->skip_merge);

  all_cmd->add_group("advanced performance tweaks", "");

  all_cmd->add_param("--minimizer-size", "size of minimizers. [4, 15]")
    ->meta("INT")
    ->def("10")
    ->checker(bc::check::f::range(4, 15))
    ->setter(options->minim_size);

  all_cmd->add_param("--minimizer-type", "minimizer type (0=lexi, 1=freq).")
    ->meta("INT")
    ->def("0")
    ->checker(bc::check::f::range(0, 1))
    ->setter(options->minim_type);

  all_cmd->add_param("--repartition-type", "minimizer repartition (0=unordered, 1=ordered).")
    ->meta("INT")
    ->def("0")
    ->checker(bc::check::f::range(0, 1))
    ->setter(options->repart_type);

  all_cmd->add_param("--nb-partitions", "number of partitions (0=auto).")
    ->meta("INT")
    ->def("0")
    ->checker(bc::check::is_number)
    ->setter(options->nb_parts);

  all_cmd->add_param("--restrict-to", "Process only a fraction of partitions. [0.05, 1.0]")
    ->meta("FLOAT")
    ->def("1.0")
    ->checker(bc::check::f::range(0.05, 1.0))
    ->setter(options->restrict_to);

  auto rtl_setter = [options](const std::string& v) {
    auto partitions = bc::utils::split(v, ',');
    for (auto& p : partitions)
      options->restrict_to_list.push_back(bc::utils::lexical_cast<uint32_t>(p));
  };

  all_cmd->add_param("--restrict-to-list", "Process only some partitions, comma separated.")
    ->meta("STR")
    ->def("")
    ->setter_c(rtl_setter);

  all_cmd->add_param("--focus", "0: focus on disk usage, 1: focus on speed. [0.0, 1.0]")
    ->meta("FLOAT")
    ->def("0.5")
    ->checker(bc::check::f::range(0.0, 1.0))
    ->setter(options->focus);

  all_cmd->add_param("--cpr", "compression for kmtricks's tmp files.")
    ->as_flag()
    ->setter(options->lz4);

  all_cmd->add_group("hash mode configuration", "");

  all_cmd->add_param("--bloom-size", "bloom filter size")
    ->meta("INT")
    ->def("10000000")
    ->checker(bc::check::is_number)
    ->setter(options->bloom_size);

  auto format_setter = [options](const std::string& v) {
    options->out_format = str_to_format(v);
  };

  all_cmd->add_param("--bf-format", "bloom filter format. [howdesbt|sdsl]")
    ->meta("STR")
    ->def("howdesbt")
    ->checker(bc::check::f::in("howdesbt|sdsl"))
    ->setter_c(format_setter);

  add_common(all_cmd, options);

  return options;
}

km_options_t repart_cli(std::shared_ptr<bc::Parser<1>> cli, repart_options_t options)
{
  bc::cmd_t repart_cmd = cli->add_command("repart", "Compute minimizer repartition.");

  repart_cmd->add_param("--file", "kmtricks input file, see README.md.")
    ->meta("FILE")
    ->checker(bc::check::is_file)
    ->setter(options->fof);

  repart_cmd->add_param("--run-dir", "kmtricks runtime directory.")
    ->meta("DIR")
    ->def("km_dir")
    ->checker(dir_already_exists)
    ->setter(options->dir);

  repart_cmd->add_param("--kmer-size", fmt::format("size of a k-mer. [8, {}]", KL[KMER_N-1]-1))
    ->meta("INT")
    ->def("31")
    ->checker(bc::check::f::range(8, KL[KMER_N-1]-1))
    ->setter(options->kmer_size);

  repart_cmd->add_group("advanced performance tweaks", "");

  repart_cmd->add_param("--minimizer-size", "size of minimizers. [4, 15]")
    ->meta("INT")
    ->def("10")
    ->checker(bc::check::f::range(4, 15))
    ->setter(options->minim_size);

  repart_cmd->add_param("--minimizer-type", "minimizer type (0=lexi, 1=freq).")
    ->meta("INT")
    ->def("0")
    ->checker(bc::check::f::range(0, 1))
    ->setter(options->minim_type);

  repart_cmd->add_param("--repartition-type", "minimizer repartition (0=unordered, 1=ordered).")
    ->meta("INT")
    ->def("0")
    ->checker(bc::check::f::range(0, 1))
    ->setter(options->repart_type);

  repart_cmd->add_param("--nb-partitions", "number of partitions (0=auto).")
    ->meta("INT")
    ->def("0")
    ->checker(bc::check::is_number)
    ->setter(options->nb_parts);

  repart_cmd->add_param("--bloom-size", "bloom filter size")
    ->meta("INT")
    ->def("10000000")
    ->checker(bc::check::is_number)
    ->setter(options->bloom_size);

  add_common(repart_cmd, options);

  return options;
}

km_options_t superk_cli(std::shared_ptr<bc::Parser<1>> cli, superk_options_t options)
{
  bc::cmd_t superk_cmd = cli->add_command("superk", "Compute super-k-mers.");

  superk_cmd->add_param("--run-dir", "kmtricks runtime directory.")
    ->meta("DIR")
    ->setter(options->dir);

  superk_cmd->add_param("--id", "sample ID, as define in the input fof.")
    ->meta("STR")
    ->setter(options->id);

  auto rtl_setter = [options](const std::string& v) {
    auto partitions = bc::utils::split(v, ',');
    for (auto& p : partitions)
      options->restrict_to_list.push_back(bc::utils::lexical_cast<uint32_t>(p));
  };

  superk_cmd->add_param("--restrict-to-list", "process only some partitions, comma separated.")
    ->meta("STR")
    ->def("")
    ->setter_c(rtl_setter);

  superk_cmd->add_param("--cpr", "output compressed super-k-mers.")
    ->as_flag()
    ->setter(options->lz4);

  add_common(superk_cmd, options);

  return options;
}

km_options_t count_cli(std::shared_ptr<bc::Parser<1>> cli, count_options_t options)
{
  bc::cmd_t count_cmd = cli->add_command("count", "Count k-mers/hashes in partitions.");

  count_cmd->add_param("--id", "sample ID, as define in kmtricks fof.")
    ->meta("STR")
    ->setter(options->id);

  count_cmd->add_param("--run-dir", "kmtricks runtime directory.")
    ->meta("DIR")
    ->checker(bc::check::is_dir)
    ->setter(options->dir);

  count_cmd->add_param("--abundance-min", "min abundance to keep a k-mer/hash.")
    ->meta("INT")
    ->def("2")
    ->checker(bc::check::is_number)
    ->setter(options->c_ab_min);

  count_cmd->add_param("--partition-id", "partition id (default: all partitions are processed.")
    ->meta("INT")
    ->def("-1")
    ->checker(bc::check::is_number)
    ->setter(options->partition_id);

  count_cmd->add_param("--mode", "count k-mers or hashes. [kmer|hash|vector|kff|kff-sk]")
    ->meta("STR")
    ->checker(bc::check::f::in("kmer|hash|vector|kff|kff-sk"))
    ->setter(options->format);

  count_cmd->add_param("--hist", "compute k-mer histograms.")
    ->as_flag()
    ->setter(options->hist);

  count_cmd->add_param("--clear", "clear super-k-mer files.")
    ->as_flag()
    ->setter(options->clear);

  count_cmd->add_param("--cpr", "output compressed partitions.")
    ->as_flag()
    ->setter(options->lz4);

  add_common(count_cmd, options);
  return options;
}

km_options_t merge_cli(std::shared_ptr<bc::Parser<1>> cli, merge_options_t options)
{
  bc::cmd_t merge_cmd = cli->add_command("merge", "Merge partitions.");

  merge_cmd->add_param("--run-dir", "kmtricks runtime directory.")
    ->meta("DIR")
    ->checker(bc::check::is_dir)
    ->setter(options->dir);

  merge_cmd->add_param("--partition-id", "partition id (-1 = all partitions are processed).")
    ->meta("INT")
    ->def("-1")
    ->checker(bc::check::is_number)
    ->setter(options->partition_id);

  auto a_min_setter = [options](const std::string& v) {
    auto [f, _] = bc::check::is_file("", v);
    if (f)
    {
      options->m_ab_min_path = v;
      return;
    }
    if (v.find('.') != std::string::npos)
    {
      bc::check::throw_if_false(bc::check::f::range(0.0, 1.0)("--abundance-min<float>", v));
      options->m_ab_min_f;
      options->m_ab_float = true;
    }
    else
      options->m_ab_min = bc::utils::lexical_cast<uint32_t>(v);
  };

  merge_cmd->add_param("--abundance-min", "min abundance to keep a k-mer/hash, see README.")
    ->meta("INT/STR/FLOAT")
    ->def("1")
    ->setter_c(a_min_setter);

  merge_cmd->add_param("--recurrence-min", "min recurrence to keep a k-mer/hash.")
    ->meta("INT")
    ->def("1")
    ->checker(bc::check::is_number)
    ->setter(options->r_min);

  merge_cmd->add_param("--save-if", "save a non-solid k-mer if it occurs in N other samples.")
    ->meta("INT")
    ->def("0")
    ->checker(bc::check::is_number)
    ->setter(options->save_if);

  auto mode_setter = [options](const std::string& v) {
    auto s = bc::utils::split(v, ':');
    options->count_format = str_to_cformat(s[0]);
    options->mode = str_to_mode(s[1]);
    options->format = str_to_format2(s[2]);
  };

  auto mode_checker = [](const std::string& p, const std::string& v) -> bc::check::checker_ret_t {
    std::string available = "kmer:pa:text|"
                            "kmer:pa:bin|"
                            "kmer:count:text|"
                            "kmer:count:bin|"
                            "hash:count:text|"
                            "hash:count:bin|"
                            "hash:pa:text|"
                            "hash:pa:bin|"
                            "hash:bf:bin|"
                            "hash:bft:bin";
    auto s = bc::utils::split(v, ':');
    std::string mode;
    std::string format;
    std::string out;

    if (s.size() != 3)
      goto fail;

    mode = s[0];
    format = s[1];
    out = s[2];

    if (format != "count" && format != "pa" && format != "bf" && format != "bft")
      goto fail;

    if (mode == "kmer")
    {
      if (format == "bf" || format == "bft")
        goto fail;
    }
    else if (mode != "hash")
      goto fail;

    if ((format == "bf" || format == "bft") && out == "text")
      goto fail;

    goto success;

    fail:
      return std::make_tuple(false, bc::utils::format_error(
        p, v, fmt::format("Available -> {}", available)));

    success:
      return std::make_tuple(true, "");
  };

  merge_cmd->add_param("--mode", "matrix mode <mode:format:out>, see README")
    ->meta("MODE:FORMAT:OUT")
    ->def("kmer:count:bin")
    ->checker(mode_checker)
    ->setter_c(mode_setter);

  merge_cmd->add_param("--clear", "clear partition files.")
    ->as_flag()
    ->setter(options->clear);

  merge_cmd->add_param("--cpr", "output compressed matrices.")
    ->as_flag()
    ->setter(options->lz4);

  add_common(merge_cmd, options);
  return options;
}

km_options_t format_cli(std::shared_ptr<bc::Parser<1>> cli, format_options_t options)
{
  bc::cmd_t format_cmd = cli->add_command("format", "Convert Bloom matrices to sample-specific Bloom filters.");

  format_cmd->add_param("--run-dir", "kmtricks runtime directory.")
    ->meta("DIR")
    ->checker(bc::check::is_dir)
    ->setter(options->dir);

  auto format_setter = [options](const std::string& v) {
    options->out_format = str_to_format(v);
  };

  format_cmd->add_param("--id", "Sample id")
    ->meta("STR")
    ->def("all")
    ->setter(options->id);

  format_cmd->add_param("--out-format", "bloom filters format. [sdsl|howdesbt]")
    ->meta("STR")
    ->def("howdesbt")
    ->checker(bc::check::f::in("howdesbt|sdsl"))
    ->setter_c(format_setter);

  format_cmd->add_param("--from-vec", "build bloom filters from bit-vectors. (kmtricks pipeline --skip-merge)")
    ->as_flag()
    ->setter(options->from_vec);

  format_cmd->add_param("--from-hash", "build bloom filters from counted hashes.")
    ->as_flag()
    ->setter(options->from_hash);

  format_cmd->add_param("--cpr-in", "compressed inputs.")
    ->as_flag()
    ->setter(options->lz4);

  format_cmd->add_param("--clear", "clear input files.")
    ->as_flag()
    ->setter(options->clear);

  add_common(format_cmd, options);
  return options;
}

km_options_t dump_cli(std::shared_ptr<bc::Parser<1>> cli, dump_options_t options)
{
  bc::cmd_t dump_cmd = cli->add_command("dump", "Dump kmtricks's files in human readable format.");

  dump_cmd->add_param("--run-dir", "kmtricks runtime directory")
    ->meta("DIR")
    ->checker(bc::check::is_dir)
    ->setter(options->dir);

  dump_cmd->add_param("--input", "path to file.")
    ->meta("FILE")
    ->checker(bc::check::is_file)
    ->setter(options->input);

  dump_cmd->add_param("-o/--output", "output file.")
    ->meta("FILE")
    ->def("stdout")
    ->setter(options->output);

  add_common(dump_cmd, options);
  return options;
}

km_options_t agg_cli(std::shared_ptr<bc::Parser<1>> cli, agg_options_t options)
{
  bc::cmd_t agg_cmd = cli->add_command("aggregate", "Aggregate partition files.");
  agg_cmd->add_param("--run-dir", "kmtricks runtime directory.")
    ->meta("DIR")
    ->checker(bc::check::is_dir)
    ->setter(options->dir);

  agg_cmd->add_group("file type", "");

  auto check_count = [](const std::string& s, const std::string& v) -> bc::check::checker_ret_t {
    if (v.find(':') == std::string::npos)
      return std::make_tuple(false, "Invalid option format.");
    std::string type = bc::utils::split(v, ':')[1];
    return bc::check::f::in("kmer|hash")(s, type);
  };

  auto set_count = [options](const std::string& v) {
    options->id = bc::utils::split(v, ':')[0];
    options->count = bc::utils::split(v, ':')[1];
  };

  agg_cmd->add_param("--count", "aggregate counted k-mers/hashes. [id:kmer|hash]")
    ->meta("ID:TYPE")
    ->def("")
    ->checker(check_count)
    ->setter_c(set_count);


  agg_cmd->add_param("--matrix", "aggregate count matrices. [kmer|hash]")
    ->meta("TYPE")
    ->def("")
    ->checker(bc::check::f::in("kmer|hash"))
    ->setter(options->matrix);


  agg_cmd->add_param("--pa-matrix", "aggregate presence/absence matrices. [kmer|hash]")
    ->meta("TYPE:P")
    ->def("")
    ->checker(bc::check::f::in("kmer|hash"))
    ->setter(options->pa_matrix);

  agg_cmd->add_group("I/O options", "");

  agg_cmd->add_param("--format", "dump in human readable format. [text|bin]")
    ->meta("STR")
    ->def("text")
    ->checker(bc::check::f::in("text|bin"))
    ->setter(options->format);

  agg_cmd->add_param("--sorted", "sorted output (A < C < T < G).")
    ->as_flag()
    ->setter(options->sorted);

  agg_cmd->add_param("--cpr-in", "compressed inputs.")
    ->as_flag()
    ->setter(options->lz4_in);

  agg_cmd->add_param("--cpr-out", "compressed output (ignored with --format text).")
    ->as_flag()
    ->setter(options->lz4);

  agg_cmd->add_param("--output", "output path.")
    ->meta("FILE")
    ->def("stdout")
    ->setter(options->output);

  add_common(agg_cmd, options);
  return options;
}

km_options_t index_cli(std::shared_ptr<bc::Parser<1>> cli, index_options_t options)
{
  bc::cmd_t index_cmd = cli->add_command("index", "Build HowDeSBT index.");
  index_cmd->add_param("--run-dir", "kmtricks runtime directory")
    ->meta("DIR")
    ->setter(options->dir);

  index_cmd->add_group("Clustering options", "");
  index_cmd->add_param("--bits", "Number of bits to use from each filter for topology computation.")
    ->meta("INT")
    ->def("100000")
    ->checker(bc::check::is_number)
    ->setter(options->bits);

  index_cmd->add_param("--cull",
                       "Remove nodes for which saturation of determined is less than cull")
    ->meta("FLOAT")
    ->def("0.2")
    ->checker(bc::check::f::range(0.0, 1.0))
    ->setter(options->cull);

  index_cmd->add_param("--keep-all-nodes", "keep all node of the binary tree.")
    ->as_flag()
    ->setter(options->keep);

  index_cmd->add_group("Build options", "");

  index_cmd->add_param("--howde", "equivalent to --determined,brief --rrr")
    ->as_flag()
    ->setter(options->howde);

  index_cmd->add_param("--allsome", "create tree nodes as all/some bloom filters")
    ->as_flag()
    ->setter(options->allsome);

  index_cmd->add_param("--determined", "create tree nodes as determined/how bloom filters")
    ->as_flag()
    ->setter(options->determined);

  index_cmd->add_param("--determined-brief",
                       "create tree nodes as determined/how, but only store  active bits")
    ->as_flag()
    ->setter(options->brief);

  index_cmd->add_param("--uncompressed", "create the nodes as uncompressed bit vector(s)")
    ->as_flag()
    ->setter(options->uncompressed);

  index_cmd->add_param("--rrr", "create the nodes as rrr-compressed bit vector(s)")
    ->as_flag()
    ->setter(options->rrr);

  index_cmd->add_param("--roar", "create the nodes as roar-compressed bit vector(s)")
    ->as_flag()
    ->setter(options->roar);

  add_common(index_cmd, options);
  return options;
}

km_options_t query_cli(std::shared_ptr<bc::Parser<1>> cli, query_options_t options)
{
  bc::cmd_t query_cmd = cli->add_command("query", "Query HowDeSBT index.");
  query_cmd->add_param("--run-dir", "kmtricks runtime directory")
    ->meta("DIR")
    ->setter(options->dir);

  query_cmd->add_param("--query", "query file (fasta/fastq/txt)")
    ->meta("FILE")
    ->checker(bc::check::is_file)
    ->setter(options->query);

  query_cmd->add_param("--output", "output file.")
    ->meta("FILE")
    ->def("stdout")
    ->setter(options->output);

  query_cmd->add_param("--threshold",
                       "fraction of query kmers that must be present in a leaf to be considered a match")
    ->meta("FLOAT")
    ->def("0.7")
    ->checker(bc::check::f::range(0.0, 1.0))
    ->setter(options->threshold);

  query_cmd->add_param("--sort", "sort matched leaves by the number of query k-mers present")
    ->as_flag()
    ->setter(options->sort);

  query_cmd->add_param("--consistency-check", "check bloom filter properties across the tree")
    ->as_flag()
    ->setter(options->check);

  add_common(query_cmd, options);
  return options;
}

void info_cli(std::shared_ptr<bc::Parser<1>> cli)
{
  bc::cmd_t info_cmd = cli->add_command("infos", "Show version and build infos.");
}

};