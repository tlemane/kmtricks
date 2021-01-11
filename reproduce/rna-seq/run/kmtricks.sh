#!/bin/bash

in_fof=$1

kmtricks.py run --file ${in_fof} --run-dir ./run --count-abundance-min 2 \
--recurrence-min 1 --mode bf_trp --split howde --max-hash 2000000000 \
--nb-cores 20 --kmer-size 20 --max-memory 1000 --nb-partitions 300 --lz4

ls *./run/storage/vectors/howde/*.bf > filterlist

howdesbt cluster --list=filterlist \
--cull=20% --bits=500K --tree=uncompressed.culled.sbt --nodename=node{number}.bf
howdesbt build --determined,brief \
--rrr --tree=uncompressed.culled.sbt --outtree=howde.culled.rrr.sbt
