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
rm -rf ./count_matrix_run
kmtricks.py run --file ./data/fof.txt \
                --run-dir ./count_matrix_run \
                --kmer-size 20 \
                --nb-cores 8 \
                --nb-partitions 4 \
                --count-abundance-min 1 \
                --recurrence-min 1 \
                --mode ascii \
                --lz4