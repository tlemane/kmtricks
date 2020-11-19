#!/bin/bash

# kmtricks pipeline: non-regression testing

function check_exit_code {
    if [ $? -ne 0 ]; then
        exit 1;
    fi
}

rm -rf ./run_test

binary=../../bin/kmtricks
fof=./data/fof.txt
nbpart=4
kmersize=20
cores=12

mkdir -p run_test/storage
cp -r ./data/partition_storage_gatb ./run_test/storage
python3 ../../kmtricks.py run --file ${fof} --run-dir ./run_test --kmer-size ${kmersize} --mode ascii --nb-cores ${cores} --nb-partitions ${nbpart} --keep-tmp --count-abundance-min 1 --recurrence-min 1 --max-memory 1000 --lz4

#check superkmers partitions
echo -ne "Superkmers partitions ..."
for ((j=1; j<3; j++))
do
    for (( i=0; i<${nbpart}; i++ ))
    do
        f1="./res/superkparts/D${j}.superk/superKparts.${i}"
        f2="./run_test/storage/superk_partitions/D${j}.superk/superKparts.${i}"
        diff ${f1} ${f2}
        check_exit_code
    done
done

echo -e "\rSuperkmers partitions ... ok"
#
##check hashes partitions
echo -ne "Hash/Kmer partitions ..."
for (( i=0; i<${nbpart}; i++ ))
do
    for ((j=1; j<3; j++)) 
    do
        diff ./res/kmerparts/partition_${i} ./run_test/storage/kmers_partitions/partition_${i}/D${j}.kmer.lz4
        check_exit_code
    done
done
echo -e "\rHash/Kmer partitions ... ok"

#check merge and trp
echo -ne "Merge ..."
for (( i=0; i<${nbpart}; i++ ))
do
    diff ./res/matrices/partition_${i}/ascii_matrix${i}.mat ./run_test/storage/matrix/partition_${i}/ascii_matrix${i}.mat
    check_exit_code
done
echo -e "\rMerge ... ok"

#rm -rf ./run_test
