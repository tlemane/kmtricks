# kmtricks

This branch is a POC for [kff k-mer compaction](https://github.com/Kmer-File-Format). Don't use it for other things.


The procedure consists of two steps:
1. Split reads into an array of super-k-mers.
2. Read each super-k-mers and associate a count to its k-mers thanks to a k-mer count table (previously computed). If a k-mer has already been seen in another super-k-mer, the current super-k-mer is splitted into super-k-mers that contain only new k-mers.

Super-k-mers are then stored in kff format according to their minimizers. In minimizer sections, each super-k-mer is represented by their 2-bit encoding without minimizer, the minimizer position and an abundance vector.


## Installation

./install.sh -c 1 -t 0 -k "32 64" -m

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

## KFF paper

### Gallus

**Accessions**
```
SRR105788
SRR105789
SRR105792
SRR105794
SRR197985
SRR197986
```

**Input: gallus.fof**
```
D1: SRR105788_1.fastq.gz ; SRR105788_2.fastq.gz ; SRR105789_1.fastq.gz ; SRR105789_2.fastq.gz ; SRR105792_1.fastq.gz ; SRR105792_2.fastq.gz ; SRR105794_1.fastq.gz ; SRR105794_2.fastq.gz ; SRR197985_1.fastq.gz ; SRR197985_2.fastq.gz ; SRR197986_1.fastq.gz ; SRR197986_2.fastq.gz
```

**CLI**
```bash
kmtricks pipeline --file gallus.fof \
                  --run-dir gallus_dir \
                  --nb-partitions 200 \
                  --count-abundance-min 1 \
                  --mode kmer:count:bin \
                  --kmer-size 32 \
                  --until count \
                  --kff-sk-output \
                  --hist
kff-tools merge -i ./gallus_dir/counts/partition_*/*.kff -o gallus_counts.kff
```

**Accessions**
### Human
```
ERR024163
ERR024164
ERR024165
ERR024166
ERR024167
ERR024168
ERR024169
ERR024170
ERR024171
ERR024172
ERR024173
ERR024174
ERR024175
ERR024176
ERR024177
ERR024178
ERR024180
ERR024183
ERR024184
ERR024185
ERR024186
```

**Input: human.fof**
```
D1 : ERR024163_1.fastq.gz ; ERR024163_2.fastq.gz ; ERR024164_1.fastq.gz ; ERR024164_2.fastq.gz ; ERR024165_1.fastq.gz ; ERR024165_2.fastq.gz ; ERR024166_1.fastq.gz ; ERR024166_2.fastq.gz ; ERR024167_1.fastq.gz ; ERR024167_2.fastq.gz ; ERR024168_1.fastq.gz ; ERR024168_2.fastq.gz ; ERR024169_1.fastq.gz ; ERR024169_2.fastq.gz ; ERR024170_1.fastq.gz ; ERR024170_2.fastq.gz ; ERR024171_1.fastq.gz ; ERR024171_2.fastq.gz ; ERR024172_1.fastq.gz ; ERR024172_2.fastq.gz ; ERR024173_1.fastq.gz ; ERR024173_2.fastq.gz ; ERR024174_1.fastq.gz ; ERR024174_2.fastq.gz ; ERR024175_1.fastq.gz ; ERR024175_2.fastq.gz ; ERR024176_1.fastq.gz ; ERR024176_2.fastq.gz ; ERR024177_1.fastq.gz ; ERR024177_2.fastq.gz ; ERR024178_1.fastq.gz ; ERR024178_2.fastq.gz ; ERR024180_1.fastq.gz ; ERR024180_2.fastq.gz ; ERR024183_1.fastq.gz ; ERR024183_2.fastq.gz ; ERR024184_1.fastq.gz ; ERR024184_2.fastq.gz ; ERR024185_1.fastq.gz ; ERR024185_2.fastq.gz ; ERR024186_1.fastq.gz ; ERR024186_2.fastq.gz
```

**CLI**
```bash
kmtricks pipeline --file human.fof \
                  --run-dir human_dir \
                  --nb-partitions 200 \
                  --count-abundance-min 1 \
                  --mode kmer:count:bin \
                  --kmer-size 32 \
                  --until count \
                  --kff-sk-output \
                  --hist
kff-tools merge -i ./human_dir/counts/partition_*/*.kff -o human_counts.kff
```
