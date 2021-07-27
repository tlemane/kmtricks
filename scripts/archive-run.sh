#!/bin/bash

function usage()
{
  echo "Usage: "
  echo "  ./archive-run.sh [-d str]"
  echo "Options: "
  echo "  -d <DIR>  -> run directory."
  echo "  -o <FILE> -> output."
}

dir=""
output=""

while getopts "d:o:" option; do
  case "$option" in
    d)
      dir=${OPTARG}
      ;;
    o)
      output=${OPTARG}
      ;;
  esac
done

tar -czvf ${output} \
          ${dir}/build_infos.txt \
          ${dir}/options.txt \
          ${dir}/hash.info \
          ${dir}/repartition_gatb \
          ${dir}/minimizers \
          ${dir}/config_gatb \
          ${dir}/superkmers/*/SuperKmerBinInfoFile \
          ${dir}/superkmers/*/PartiInfoFile \
          ${dir}/partition_infos \
          ${dir}/histograms \
          ${dir}/kmtricks.fof \
          ${dir}/merge_infos

