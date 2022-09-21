# kmtricks

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

![kmtricks](https://github.com/tlemane/kmtricks/workflows/kmtricks/badge.svg)

kmtricks is a modular tool suite for counting kmers, and constructing Bloom filters or kmer matrices, for large collections of sequencing data.

## Rationale

kmtricks is optimized for the analysis of **multiple FASTA/FASTQ** files (gzipped or not). It features:
 * Fast **k-mer matrix** construction
 * Fast **Bloom filters** construction
 * **Rescues low-abundance k-mers** when they are seen in multiple samples

Note: for counting kmers from a single file, kmtricks works but is slightly slower than a traditional k-mer counter (e.g. KMC). It is really optimized for merging count information across multiple samples, which traditional k-mer counters cannot do.

## Overview

**Input**: a set of read sets in FASTA or FASTQ format, gzipped or not.

**Final output** is either:

* a matrix of kmer abundances. M<sub>i,j</sub> is the abundance of kmer i in the read set j
* a matrix of kmer membership. M<sub>i,j</sub> is the presence (1) or absence (0) of kmer i in the read set j
* a vector of Bloom filters. M<sub>i,j</sub> is the presence (1) or absence (0) of the hash_value i (line numbers are hash values) in the read set j.
  * In this case, this matrix is provided vertically (one column is a bloom filter corresponding to one dataset).
  * After transposition, this matrix may also be provided horizontally (one line is a Bloom filter corresponding to one dataset). This enables to provide efficiently an independent Bloom filter per input read file.

## Installation and usage

Instructions for installation and usage are provided in the [wiki](https://github.com/tlemane/kmtricks/wiki/Home).

## Limitations

kmtricks needs disk space to run. The disk usage is variable and depends on data, parameters and output format. Based on our observations, the required space is between 20% of the total input size (gzipped) and the total input size (including outputs).

## Reporting an issue

If you encounter a problem, please open an issue with a description of your run and the return of [kmtricks infos](https://github.com/tlemane/kmtricks/wiki/infos). If you encounter a critical error like a segmentation fault, kmtricks automatically dumps a file `kmtricks_backtrace.log` in your current directory. This file is somewhat illegible in release mode. If you can, compile kmtricks in debug mode, launch it again and join the content of this file. If you cannot directly compile kmtricks on your system, the conda package provides `kmtricks-debug` binary for this case.

## Reference


T. Lemane, P. Medvedev, R. Chikhi and P. Peterlongo, "kmtricks: Efficient and flexible construction of Bloom filters for large sequencing data collections." Bioinformatics Advances, 2022, doi:10.1093/bioadv/vbac029.
```tex
@article{kmtricks,
    author = {Lemane, Téo and Medvedev, Paul and Chikhi, Rayan and Peterlongo, Pierre},
    title = "{kmtricks: Efficient and flexible construction of Bloom filters for large sequencing data collections}",
    journal = {Bioinformatics Advances},
    year = {2022},
    doi = {10.1093/bioadv/vbac029},
    url = {https://doi.org/10.1093/bioadv/vbac029},
}
```

## Contacts

Téo Lemane: teo[dot]lemane[at]inria[dot]fr

Rayan Chikhi: rayan[dot]chikhi[at]pasteur[dot]fr

Pierre Peterlongo: pierre[dot]peterlongo[at]inria[dot]fr
