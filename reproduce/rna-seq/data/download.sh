#!/bin/bash

mkdir -p ./fasta

c=0
m=100

if [ ! -z "$1" ]; then
  m=$1
fi

while read -r exp;
do
  if [[ ! -f "./fastq/${exp}.fasta.gz" ]]; then
    parallel-fastq-dump --sra-id ${exp} --threads 8 --outdir ./fasta --gzip --fasta
  fi
  let "c++"
  if (( $c == $m )); then
    break
  fi
done < sra_id.txt
