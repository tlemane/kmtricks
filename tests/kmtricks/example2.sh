#!/bin/bash

check_bin(){
    which $1 &> /dev/null
    if [ $? -ne 0 ]
    then
        echo "$0 requires $1$reset"
        exit 1
    fi
}

check_bin kmtricks.py

echo "> BLOOM FILTERS CONSTRUCTION <"
rm -rf ./km_index
kmtricks.py run --file ./data/fof.txt \
                --run-dir ./km_index \
                --kmer-size 20 \
                --nb-cores 8 \
                --nb-partitions 4 \
                --count-abundance-min 1 \
                --recurrence-min 1 \
                --mode bf_trp \
                --hasher sabuhash \
                --max-hash 100000 \
                --split howde \
                --lz4

echo "\n> INDEX CONSTRUCTION <"
cd ./km_index/storage/vectors/howde
ls *.bf > bf_list
km_howdesbt cluster bf_list
km_howdesbt build --howde bf_list

echo "\n> QUERIES <"
cp ../../../../data/1.fasta queries.fa
echo ">RandomSeq" >> ./queries.fa
echo "AGCAGCACTACTACACATCATCTATTCACACAAGCGGGATATAGTAGAGATATAAGCAGGCAGACGACGCAGCAGC" >> ./queries.fa
km_howdesbt queryKm \
            --tree=bf_list.detbrief.sbt \
            --repart=../../partition_storage_gatb/minimRepart.minimRepart \
            --win=../../hash_window.vec \
            queries.fa