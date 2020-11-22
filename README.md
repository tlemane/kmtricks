# kmtricks

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

| linux            | osx             |
|------------------|-----------------|
| [![Linux][1]][3] | [![osx][2]][3]  |

[1]: https://travis-matrix-badges.herokuapp.com/repos/tlemane/kmtricks/branches/master/1?
[2]: https://travis-matrix-badges.herokuapp.com/repos/tlemane/kmtricks/branches/master/2?
[3]: https://travis-ci.org/github/tlemane/kmtricks



kmtricks is a tool suite for counting kmers, and constructing bloom filters or counted kmer matrices from large and numerous read sets. 

## kmtricks IOs

kmtricks is composed of a **set of independent modules** designed for kmer counting given a set of raw read sets. 

A pipeline of those modules is proposed, with the following IOs. See the [kmtricks pipeline](#kmtricks-pipeline) section for details.

**Input** is composed of a set of read sets in fasta or fastq format, gzipped or not.

**Final output** depends on the user usage. It may be

* a matrix of kmer x abundance. M_(i,j) is the abundance of kmer i in the read set j
* a matrix of kmer x presence or absence. M_(i,j) is the presence (1) or absence (0) of kmer i in the read set j
* a matrix of bloom filters. M_(i,j) is the presence (1) or absence (0) of the hash_value i (line numbers are hash values) in the read set j.
  * In this case, this matrix is provided vertically (one column is a bloom filter corresponding to one dataset).
  * After transposition, this matrix may also be provided horizontally (one line is a bloom filter corresponding to one dataset). This enables to provide efficiently an independent bloom filter per input read file.  

## kmtricks performances

Compared to a usual pipeline as the one used by `HowDeSBT` using `JellyFish` for generating a bloom filter per input read set, `kmtricks` is 4.2 times faster.

|      Method     | Indexation time | Max memory | Disk usage |
|-----------------|:---------------:|:----------:|:----------:|
| HowDeSBT makebf |       2h27      |   13.2 GB  |   55.1 GB  |
| kmtricks        |      35min48s   |   3.5 GB   |   56.6 GB  

Test realised with 100 RNA-seq experiments with 20 cores 100 GB RAM Intel(R) Xeon(R) 2.60GHz

List of IDs available [here](tests/kmtricks/experiment_list_100.txt).

## kmtricks usage

kmtricks can be used in two different ways: by using each **independent modules** or by using the **pipeline** (kmtricks binary). 

### kmtricks modules

kmtricks is composed of 5 independent modules

<img src="https://github.com/tlemane/kmtricks/blob/master/doc/kmtricks_pipeline.png" width="500">  

**Note1:** Using any of those modules requires the existence of the `run-dir` directory and its whole internal structure. The creation of the directory and its structure can be done thanks to the following command: 

`python3 kmtricks.py env`

```
my_run_directory/  
 ├── logs  
 │   ├── cmds.log  
 │   ├── counter  
 │   ├── merger  
 │   ├── superk  
 │   └── split  
 └── storage
     ├── config_storage_gatb // gatb config
     │   └── config.config
     ├── fof.txt        
     ├── hash_window.vec     
     ├── kmers_partitions    // km_superk_to_kmer_count output
     │   ├── partition_0
     │   └── partition_1
     ├── matrix              // km_merge_within_partition output
     │   ├── partition_0
     │   └── partition_1
     ├── superk_partitions   // km_reads_to_superk output
     └── vectors             // km_output_convert output 
         ├── howde
         └── sdsl
```

**Note2:** Run any of the binary with no argument provides a detailed list of options and mandatory arguments.

Each module is presented below. However, the `kmtricks` binary enables to execute automatically all modules. See the [kmtricks pipeline](#kmtricks-pipeline) section.

#### Module `km_minim_repart`: determine partitions

From reads, determine minimizers and assign each minimizer to a partition.

Example: `./bin/km_minim_repart -file file_of_files.txt -kmer-size 31 -run-dir my_directory_output_name`

#### Module `km_reads_to_superk`: from reads to partitioned super kmers

For each read file,  using the previously determined partitions from minimizers, write superkmers into corresponding partitions

Example: `./bin/km_reads_to_superk -file read_file.fasta -run-dir my_directory_output_name -nb-cores 8 -kmer-size 31`

#### Module `km_superk_to_kmer_counts`: from super kmers to counted elements 

For one superkmer partition, determine, sort and count elements that may be kmers or hash value.

Example: `./bin/km_reads_to_kmer_counts -file read_file.fasta -run-dir my_directory_output_name -kmer-size 31 -part-id N`

Option `-mode` enables to provide results either as kmers or hash values 

#### Module `km_merge_within_partition ` merges counted kmers and transpose matrix

For a given partition id, merges values for all input read files. 

Example: `./bin/km_merge_within_partition -run-dir my_directory_output_name -part-id 0 -abundance-min 2 -recurrence-min 2`

#### Module `km_output_convert`: generates output for downstream usages

Given the merged partitions, depending on the user choice, outputs a SDSL compatible or a HowDeSBT compatible set of files. 

Example: `./bin/km_output_convert -run-dir my_directory_output_name -nb-files nb_of_reads_files -split howde -kmer-size 31`

### kmtricks pipeline

`kmtricks.py` is a pipeline of the five modules.


Note that this binary also enables to run independently any module (option `--only`) or enables to run modules until a step (option `--until`).

**Usage:**

`python3 kmtricks.py run --file file_of_files.txt --run-dir my_directory_output_name`

Final results are stored in the `directory_output_name/storage/matrix/` and `directory_output_name/storage/vectors/`.

```
usage: kmtricks.py run --file FILE --run-dir DIR [--kmer-size INT] [--count-abundance-min INT]
                       [--abundance-max INT] [--max-count INT] [--max-memory INT] [--mode STR]
                       [--nb-cores INT] [--merge-abundance-min INT/STR] [--recurrence-min INT]
                       [--save-if INT] [--skip-merge] [--until STR] [--only STR]
                       [--minimizer-type INT] [--minimizer-size INT] [--repartition-type INT]
                       [--nb-partitions INT] [--hasher STR] [--max-hash INT] [--split STR]
                       [--keep-tmp] [--lz4] [-h]

kmtricks pipeline

global:
  --file FILE                    fof that contains path of read files, one per line [required]
  --run-dir DIR                  directory to write tmp and output files [required]
  --kmer-size INT                size of a kmer [default: 31]
  --count-abundance-min INT      min abundance threshold for solid kmers [default: 2]
  --abundance-max INT            max abundance threshold for solid kmers [default: 3000000000]
  --max-count INT                allows to deduce the integer size to store counts [default: 255]
  --max-memory INT               max memory available in megabytes [default: 8000]
  --mode STR                     output matrix format: [bin|ascii|pa|bf|bf_trp] [default: bin]
  --nb-cores INT                 number of cores [default: 8]
  --keep-tmp                     keep all tmp files [no arg]
  --lz4                          lz4 compression for tmp files [no arg]
  -h, --help                     Show this message and exit

merge options:
  --merge-abundance-min INT/STR  min abundance threshold for solid kmers [default: 1]
  --recurrence-min INT           min reccurence threshold for solid kmers [default: 1]
  --save-if INT                  keep a non-solid kmer if it's solid in X dataset [default: 0]
  --skip-merge                   skip merge step, only with --mode bf [no arg]

pipeline control:
  --until STR                    run until step: [repart|superk|count|merge|split] [default: all]
  --only STR                     run until step: [repart|superk|count|merge|split] [default: all]

advanced performance tweaks:
  --minimizer-type INT           minimizer type (0=lexi, 1=freq) [default: 0]
  --minimizer-size INT           size of minimizer [default: 10]
  --repartition-type INT         minimizer repartition (0=unordered, 1=ordered) [default: 0]
  --nb-partitions INT            number of partitions (0=auto) [default: 0]

hash mode configuration:
  --hasher STR                   hash function: sabuhash, xor [default: xor]
  --max-hash INT                 max hash value (0 < hash < max(int64)) [default: 1000000000]
  --split STR                    split matrix in indidual files: [sdsl|howde] (only with -mf, --mode
                                 bf_trp) [default: none]
```

If you need specific thresholds for some or all datasets, you can add these thresholds in the fof. For datasets without a specific threshold, the --count-abundance-min is used.

**File of file format:**
* STR **:** PATH_1 **;** ... **;** PATH_N **!** INT
* Dataset ID : one or more fastx.gz files (seperated by **;**) ! An optional min abundance threshold

**Fof example**:
```
A1 : /path/to/fastq_A1_1
B1 : /path/to/fastq_B1_1 ; /with/mutiple/fasta_B1_2
```

**Full example, with HowDeSBT compatibility TODO update**

```bash
ls myreads1.fq.gz myreads2.fq.gz myreads3.fg.gz > file_of_files.txt # creates the file of files
python3 kmtricks.py run --file file_of_files.txt --run-dir my_directory_output_name --mode bf_trp --hasher sabuhash --split howde
```

**logs**

All execution logs are stored in the `my_directory_output_name/logs` directory.

## kmtrick librairies

In addition to modules, the `libs/kmtricks` directory contains headers. They provide a framework for creating independent tools. The `libs/snippets` directory provides usage examples of those librairies.

## Install

Maximal size of k-mers and maximal stored counts must be set at compile time for some kmtricks binaries. 

**Available values:**
* KMER_NB_BIT=8;16;32;64;128 -> respectively for k less or equal to : 4, 8, 16, 32, 64
* COUNT_NB_BIT=8;16;32 -> respectively for max counts: 255, 65535, 4294967295.


```bash
git clone --recursive https://github.com/tlemane/kmtricks
cd kmtricks
mkdir build ; cd build

# Several example, use only one
cmake .. # Default, here KMER_NB_BIT="32;64" and COUNT_NB_BIT="8;16;32"
cmake .. -DKMER_NB_BIT="64;128" -DCOUNT_NB_BIT="32" # Select values
cmake .. -DKMER_NB_BIT=ALL -DCOUNT_NB_BIT=ALL # All available values

make -j8
```

`kmtricks.py` pipeline automatically selects the binaries to be used according to parameters or provides compilation instructions if the required binaries are missing.

## Test

```bash
cd build
cmake .. -DTEST=1
ctest CTestTestfile.cmake
```

**Warning**: kmtricks is under active development. Its results, features and performances are likely to be improved quickly over time.

## Contacts

Téo Lemane: teo.lemane@inria.fr

Rayan Chikhi: rayan.chikhi@pasteur.fr

Pierre Peterlongo: pierre.peterlongo@inria.fr

