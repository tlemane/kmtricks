# kmtricks

| linux            | osx             |
|------------------|-----------------|
| [![Linux][1]][3] | [![osx][2]][3]  | 

[1]: https://travis-matrix-badges.herokuapp.com/repos/tlemane/kmtricks/branches/master/1?use_travis_com=true
[2]: https://travis-matrix-badges.herokuapp.com/repos/tlemane/kmtricks/branches/master/2?use_travis_com=true
[3]: https://travis-ci.com/github/tlemane/kmtricks

**Warning**: kmtricks is under active development. Its results, features and performances are likely to be improved quickly over time.

## kmtricks IOs

kmtricks is composed of a **set of independent modules** designed for kmer counting given a set of raw read sets. A pipeline of those modules is proposed. 

**Input** is composed of a set of read sets in fasta or fastq format, gzipped or not.

**Final output** depends on the user usage. It may be

* a matrix of kmer x abundance. M_(i,j) is the abundance of kmer i in the read set j
* a matrix of kmer x presence or absence. M_(i,j) is the presence (1) or absence (0) of kmer i in the read set j
* a matrix of bloom filters. M_(i,j) is the presence (1) or absence (0) of the hash_value i (line numbers are hash values) in the read set j.
  * In this case, this matrix is provided vertically (one column is a bloom filter corresponding to one dataset).
  * After transposition, this matrix may also be provided horizontally (one line is a bloom filter corresponding to one dataset). This enables to provide efficiently an independent bloom filter per input read file.  

## kmtricks performances

Compared to a usual pipeline as the one used by `HowDeSBT` using `JellyFish` for generating a bloom filter per input read set, `kmtricks` is 2.1 times faster.

<img src="https://github.com/tlemane/kmtricks/blob/master/doc/perf_kmtricks_100.png" width="300">

Test realised with 100 RNA-seq experiments with 20 cores 100 GB RAM Intel(R) Xeon(R) 2.60GHz

List of IDs available [here](tests/kmtricks/experiment_list_100.txt).

## kmtricks modules

kmtricks is composed of 5 independent modules

<img src="https://github.com/tlemane/kmtricks/blob/master/doc/kmtricks_pipeline.png" width="500">

**Note1:** Run any of the binary with no argument provides a detailed list of options and mandatory arguments.

**Note2:** Using any of those modules requires the existence of the `run-dir` directory and its whole internal structure. The creation of the directory and its structure can be done thanks to the following command: 

`./bin/kmtricks env`

Each module is presented below. However, the `kmtricks` binary enables to execute automatically all modules. See the [kmtricks pipeline](#kmtricks-pipeline) section.

### Module `km_minim_repart`: determine partitions

From reads, determine minimizers and assign each minimizer to a partition.

**Usage example**

`./bin/km_minim_repart -file file_of_files.txt -kmer-size 31 -run-dir my_directory_output_name`

### Module `km_reads_to_superk`: from reads to partitioned super kmers

For each read file,  using the previously determined partitions from minimizers, write superkmers into corresponding partitions

**Usage example**

`./bin/km_reads_to_superk -file read_file.fasta -run-dir my_directory_output_name -nb-cores 8 -kmer-size 31`

### Module `km_superk_to_kmer_count`: from super kmers to counted elements 

For one superkmer partition, determine, sort and count elements that may be kmers or hash value.

`./bin/km_reads_to_kmer_count -file read_file.fasta -run-dir my_directory_output_name -kmer-size 31 -part-id N`

Option `-mode` enables to provide results either as kmers or hash values 

### Module `km_merge_within_partition ` merges counted kmers and transpose matrix

For a given partition id, merges values for all input read files. 

`./bin/km_merge_within_partition -run-dir my_directory_output_name -part-id 0 -abundance-min 2 -recurrence-min 2`

### Module `km_output_convert`: generates output for downstream usages

Given the merged partitions, depending on the user choice, outputs a SDSL compatible or a HowDeSBT compatible set of files. 

`./bin/km_output_convert -run-dir my_directory_output_name -nb-files nb_of_reads_files -split howde -kmer-size 31`

## kmtricks pipeline

The `kmtricks` executable (in the `bin` directory) is a pipeline of the four modules. 

Note that this binary also enables to run independently any module (option `-only`) or enables to run modules until a step (option `-until`).

**Usage:**

`./bin/kmtricks -file file_of_files.txt -run-dir my_directory_output_name`

Final results are stored in the `directory_output_name/storage/matrix/`

**Main options**

* kmtricks options
  * -file:     					 fof that contains path of read files, one per line
  * -run-dir:                  directory to write tmp and output files
  * -kmer-size:              size of a kmer  [default '31']
  * -abundance-min:   min abundance threshold for solid kmers  [default '2']
  * -abundance-max:  max abundance threshold for solid kmers  [default '3000000']
  * -recurrence-min:   min recurrence threshold through datasets for solid kmers  [default '2'].
    * TODO details
  * -max-memory:       max memory available in megabytes  [default '8000']
  * -matrix-fmt:            output matrix format: ascii, bin, pa, bf, bf_trp  [default 'bin']. 
    * TODO details
  * -nb-cores:               number of cores  [default '8']
* kmtricks pipeline control options
  * -until:                      run until step : part, superk, count, merge  [default 'all']
  * -only:                      run only step : part, superk, count, merge  [default 'all']
* advanced performance tweaks options
  * -minimizer-type:   minimizer type (0=lexi, 1=freq)  [default '0']
  * -minimizer-size:    size of a minimizer  [default '10']
  * -repartition-type:  minimizer repartition (0=unordered, 1=ordered)  [default '0']
    * TODO details
  * -nb-parts:               number of partitions  [default '0']
* hash mode configuration, only with -matrix-fmt <bf | bf_trp> options
  * -hasher:                  hash function: sabuhash, xor  [default 'xor']. 
    * For compatibility with HowDeSBT, use "sabuhash".
  * -max-hash:             max hash value ( 0 < hash < max(int64) )  [default '1000000000']. 
    * This is also the size of the final bloom filters 
  * -split:                       split matrix in individual files: sdsl, howde, (only with -matrix-fmt bf_trp)  [default 'none']. 
    * sdsl: dump one sdsl::bit_vector per file
    * howde: dump one bloom filter per file, with HowDeSBT compatibility

**Full example, with HowDeSBT compatibility**

```bash
ls myreads1.fq.gz myreads2.fq.gz myreads3.fg.gz > file_of_files.txt # creates the file of files
./bin/kmtricks -file file_of_files.txt -run-dir my_directory_output_name -matrix-fmt bf_trp -hasher sabuhash -split howde
```

**logs**

All execution logs are stored in the `my_directory_output_name/logs` directory.

## kmtrick librairies

In addition to modules, the `libs/kmtricks` directory contains headers `merger.hpp` and `bitmatrix.hpp`. They provide a framework for creating independent tools. The `libs/snippets` directory provides usage examples of those librairies.

## Install

```bash
git clone --recursive https://github.com/tlemane/kmtricks
cd kmtricks
sh install.sh
```

## Test

```bash
cd build
ctest CTestTestfile.cmake
```

## Contacts

Téo Lemane: teo.lemane@inria.fr

Rayan Chikhi: rayan.chikhi@univ-lille.fr

Pierre Peterlongo: pierre.peterlongo@inria.fr

