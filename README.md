# kmtricks
[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

|                                                              linux                                                             |                                                               osx                                                              |
|:------------------------------------------------------------------------------------------------------------------------------:|:------------------------------------------------------------------------------------------------------------------------------:|
| [![Travis](https://api.travis-ci.com/tlemane/kmtricks.svg?branch=master&job=1)](https://travis-ci.com/github/tlemane/kmtricks) | [![Travis](https://api.travis-ci.com/tlemane/kmtricks.svg?branch=master&job=4)](https://travis-ci.com/github/tlemane/kmtricks) |

**Warning**: kmtricks is under active development. Its results, features and performances are likely to be improved quickly over time.

## kmtrics IOs

kmtricks is composed of a **set of independent modules** designed for kmer counting given a set of raw read sets. A pipeline of those modules is proposed. 

**Input** is composed of a set of read sets in fasta or fastq format, gzipped or not.

**Final output** depends on the user usage. It may be

* a matrix of kmer x abundance. M_(i,j) is the abundance of kmer i in the read set j
* a matrix of kmer x presence or absence. M_(i,j) is the presence (1) or absence (0) of kmer i in the read set j
* a matrix of hash_value x abundance. M_(i,j) is the abundance of the hash_value i (line numbers are hash values) in the read set j.
* a matrix of bloom filters. M_(i,j) is the presence (1) or absence (0) of the hash_value i (line numbers are hash values) in the read set j.

## kmtrics performances

Compared to a usual pipeline as the one used by `HowDeSBT` using `JellyFish` for generating a bloom filter per input read set, `kmtricks` is 2.1 times faster.

<img src="https://github.com/tlemane/kmtricks/blob/master/doc/perf_kmtricks_100.png" width="300">

Test realised with 100 RNA-seq experiments with 20 cores 100 GB RAM Intel(R) Xeon(R) 2.60GHz

List of IDs available [here](tests/kmtricks/experiment_list_100.txt).

## kmtrics modules

kmtricks is composed of 5 independent modules

<img src="https://github.com/tlemane/kmtricks/blob/master/doc/kmtricks_pipeline.png" width="500">

### Determine partitions: `km_minim_repart`

From reads, determine minimizers and assign each minimizer to a partition

### Reads to partitioned super kmers: `km_reads_to_superk`

For each read file,  using the previously determined partitions from minimizers, write superkmers into corresponding partitions

### Super kmers to counted elements: `km_superk_to_kmer_count` 

For each superkmer partition, determine, sort and count kmers, or hash value.

(element = kmer or hash value)

### Merging counted kmers and transpose matrix: `km_merge_within_partition`

For a given partition id, merges values for all input read files. 

### Generate output for downstream usages: `km_output_convert`

Given the merged partitions, depending on the user choice, outputs a SDSL compatible or a HowDeSBT compatible set of files. 

## kmtrics pipeline

The kmtricks executable (in the `bin` directory) is a pipeline of the four modules. 

## Install

```bash
git clone --recursive https://github.com/tlemane/kmtricks
cd kmtricks
sh install.sh
```

## Test

```bash
cd build
ctest --verbose CTestTestfile.cmake
```

## Contacts

Téo Lemane: teo.lemane@inria.fr

