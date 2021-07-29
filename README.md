# kmtricks

This branch is a POC for [kff k-mer compaction](). Don't use it for other things.


The procedure consists of two steps:
1. Split reads into an array of super-k-mers.
2. Read each super-k-mers and associate a count to its k-mers thanks to a k-mer count table (previously computed). If a k-mer has already been seen in another super-k-mer, the current super-k-mer is splitted into super-k-mers that contain only new k-mers.

Super-k-mers are then stored in kff format according to their minimizers. In minimizer sections, each super-k-mer is represented by their 2-bit encoding without minimizer, the minimizer position and an abundance vector.


## Installation

See [kmtricks installation]().

## Usage

With `fof.txt`:
```
sample1: SRR_X_1.fastq.gz ; SRR_X_2.fastq.gz ; SRR_Y_1.fastq.gz ; SRR_Y_2.fastq.gz
sample1: SRR_Z_1.fastq.gz ; SRR_Z_2.fastq.gz
```

```bash
kmtricks/bin/kmtricks pipeline --file fof.txt \
                               --run-dir kffsk \
                               --count-abundance-min 1 \
                               --mode kmer:count:bin \
                               --kmer-size 32 \
                               --nb-partitions \
                               --verbose info \
                               --kff-sk-output \
                               --until count \
                               --threads 20

kff-tools merge -i ./kffsk/counts/partition_*/sample1.kff -o sample1_count.kff
kff-tools merge -i ./kffsk/counts/partition_*/sample2.kff -o sample2_count.kff
```



