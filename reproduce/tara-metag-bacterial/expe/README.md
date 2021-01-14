# kmtricks on all bacterial metaG Tara



Cf [README](../data/README.md) file in the data directory for information about data. 

**Reminder:** 

- there are approx 266 billion distinct kmers in this data set. 

- The total amount of fastq.gz is 6.1TB

## Running kmtricks

kmtricks commit [b77c45d](https://github.com/tlemane/kmtricks/commit/b77c45df3a75482aad9e9a6293d838a5adea67d9)

(this was dev branch, corresponding to release 0.0.1, and enabling to select the wanted log files)

**Command line**

Here we conserved all kmers from all dataset, before to filter the low abundant one (with thesholds as defined in file `upper_rare_kmer_thresholds_10_percent.txt`) that are not rescued.

We generated bloom filters composed of 4 billion bits each. 

```bash
python3 kmtricks.py run \
  --file bact_metaG_factorized.list \ # A corresponds to the bact_metag_tara.txt file as provided in data
  --kmer-size 20 \
  --run-dir kmtricks0_0_1_metag_bact_tara \
  --count-abundance-min 1 \
  --max-count 256 \
  --max-memory 8000 \
  --mode bf_trp \
  --nb-cores 60 \
  --lz4 \
  --merge-abundance-min upper_rare_kmer_thresholds_10_percent.txt \
  --recurrence-min 1 \
  --save-if 1 \
  --log-files merge \
  --max-hash 40000000000 \
  --hasher sabuhash \
  --split howde
```



**Partitions**

- Generates 2544 partitions
- Reminder: 241 read sets

**Computation time**

1 day, 19h, 51minutes. 

**Ressources:** 

disk used 6.2 TB - (2.93 TB if kmer seen once are removed `--count-abundance-min 2`)

max memory  50 GB

## Creation of the howDeSbt index, with option --cull

Determine the tree topology (less than a minute)

```bash 
ls kmtricks0_0_1_metag_bact_tara/storage/vectors/howde/*.bf >  list_bf.txt
kmtricks/bin/howdesbt cluster list_bf.txt --nocull --out=cluster_tara 10000000..20000000
```

`
Note that we do not compute the tree topology using the first bits of the bloom filters. This is because, these first bits (first 8 million in this experiment) correspond to the first partition, which contains low complexity kmers. 

`kmtricks/bin/howdesbt build --howde cluster_tara`

- Computation Time: RUNNING
- Max Disk Used: RUNNING
- Max Mem Used: RUNNING
- Size Directory: RUNNING
- Directory contains RUNNING bf files and one sbt file



## Query of the so-created index, with 1000 reads

TODO RESTARTS

```bash
zcat /ccc/store/cont007/fg0001/fg0001/rawdata/projet_APY/AAAI/RunsSolexa/110712_BISMUTH_63A2BAAXX/APY_AAAIOSF_1_1_63A2BAAXX_clean.fastq.gz | head -n 4000 | ./fastq2fasta.py > random1000.fa
```



```bash 
~/kmtricks_3a3426f/bin/howdesbt queryKm --tree=cluster_tara.detbrief.sbt --repart=minimRepart.minimRepart --win=hash_window.vec random1000.fa
```

Time: 

real	todo
user	todo
sys	todo

[Result file](https://gitlab.inria.fr/ppeterlo/kmtricks_tara/-/blob/master/expe/TODO)

Most found read sets: 

TODO

