# kmtricks
[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

[![Travis](https://api.travis-ci.com/tlemane/kmtricks.svg?branch=master)](https://travis-ci.com/github/tlemane/kmtricks)

## kmtrics IOs

kmtricks is a **set of independent modules** designed for kmer counting given a set of raw read sets.

**Input** is composed of a set of read sets in fasta or fastq format, gzipped or not.

**Final output** depends on the user usage. It may be

* a matrix of kmer x abundance. M_(i,j) is the abundance of kmer i in the read set j
* a matrix of kmer x presence or absence. M_(i,j) is the presence (1) or absence (0) of kmer i in the read set j
* a matrix of hash_value x abundance. M_(i,j) is the abundance of the hash_value i in the read set j.
* a matrix of hash_value x presence or absence. M_(i,j) is the presence (1) or absence (0) of the hash_value i in the read set j.

## kmtrics performances

Compared to a usual pipeline as the one used by `HowDeSBT` using `JellyFish` for generating a bloom filter per input read set, `kmtricks` is 1.6 times faster.

<img src="https://github.com/tlemane/kmtricks/blob/master/doc/perf_kmtricks_100.png" width="300">

Test realised with 100 RNA-seq experiments with 20 cores 100 GB RAM Intel(R) Xeon(R) 2.60GHz

List of IDs available [here](tests/kmtricks/experiment_list_100.txt).



## kmtrics modules

kmtricks is composed of 5 independent modules

<img src="https://github.com/tlemane/kmtricks/blob/master/doc/kmtricks_pipeline.png" width="500">

### Determine partitions

From reads, determine minimizers and assign each minimizer to a partition

### Reads to partitioned super kmers

For each read file,  using the previously determined partitions from minimizers, write superkmers into corresponding partitions

### Super kmers to counted elements 

For each superkmer partition, determine, sort and count kmers, or hash value.

(element = kmer or hash value)

### Merging counted kmers 

TODO

### Transpose matrix

TODO

## INSTALL

TODO

## Tests

TODO
