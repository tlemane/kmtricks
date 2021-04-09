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
check_bin kmtricks-socks.py

kmtricks-socks.py data/fof_socks.txt -m 1
kmtricks-socks.py lookup-kmer ./data/socks_query