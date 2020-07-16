#!/bin/bash

# kmtricks pipeline: non-regression testing

function check_exit_code {
    if [ $? -ne 0 ]; then
        exit 1;
    fi
}

binary=../../bin/kmtricks
fof=./data/fof.txt
nbpart=4
kmersize=20

${binary} -file ${fof} -run-dir ./run_test -kmer-size ${kmersize} -matrix-fmt 4 -nb-cores 1 -nb-parts ${nbpart} -max-memory 1000 -keep-tmp 1 -max-hash 100 -abundance-min 1 -recurrence-min 1 #-only 4

#check superkmers partitions
echo -ne "Superkmers partitions ..."
while IFS= read -r file
do
    for (( i=0; i<${nbpart}; i++ ))
    do
        file="$(basename -- $file)"
        diff ./res/superkparts/${file}.superk/superKparts.${i} ./run_test/storage/superk_partitions/${file}.superk/superKparts.${i}
        check_exit_code
    done
done < ${fof}
echo -e "\rSuperkmers partitions ... ok"

#check hashes partitions
echo -ne "Hash/Kmer partitions ..."
for (( i=0; i<${nbpart}; i++ ))
do
    while IFS= read -r file
    do
        file="$(basename -- $file)"
        diff ./res/kmerparts/partition_${i} ./run_test/storage/kmers_partitions/partition_${i}/${file}.kmer
        check_exit_code
    done < ${fof}
done
echo -e "\rHash/Kmer partitions ... ok"

#check merge and trp
echo -ne "Merge/Tranpose ..."
for (( i=0; i<${nbpart}; i++ ))
do
    diff ./res/matrices/partition_${i}/no_trp_bf${i}.mat ./run_test/storage/matrix/partition_${i}/no_trp_bf${i}.mat
    check_exit_code
    diff ./res/matrices/partition_${i}/trp_bf${i}.mat ./run_test/storage/matrix/partition_${i}/trp_bf${i}.mat
    check_exit_code
done
echo -e "\rMerge/Transpose ... ok"

rm -rf ./run_test
