#!/bin/bash

in_fof=$1

while read -r exp;
do
  IFS=" " read name cutoff <<< "$exp"
  nf="$(basename -- $name)"
  mccortex31 build -m 40G -k 31 -n 2G -t 20 -s $nf -1 $name output/$nf.ctx &> mccortex.log
done < ${in_fof}

cobs compact-construct -T 20 -k 31 ./output index.cobs_compact
