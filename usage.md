# kmtricks usage

**About input:**  
One sample per line, with an ID, a list of files and an optional solid threshold.
`fof.txt` with count-abundance-min specified for each sample, otherwise (`D: sample.fastq.gz`) use --count-abundance-min INT from the command line (default is 2).
```
D1: sample1.fastq.gz ! 2
D2: sample2.fastq.gz ! 3
```

**About outputs:**  

```
my_run_directory/  
├── config.log # log about general parameter, git sha1, ect...
├── logs #enable all log with --debug or --log-files [repart,superk,count,merge,split]
│   ├── cmds.log # summary of module calls
│   ├── counter
│   ├── merger
│   ├── split
│   └── superk
└── storage
    ├── config_storage_gatb         # (0) GATB configuration files (kmtricks.py env)
    ├── partition_storage_gatb      # (1) --until repart (queryable using repartition.hpp)
    │   └── minimRepart.minimRepart # minimizers repartition (required for queries)
    ├── superk_partitions           # (2) --until superk (readable using skreader.hpp)
    ├── kmers_partitions            # (3) --until count  (stream matrix from here using merger.hpp)
    ├── matrix                      # (4) --until merge --mode [ascii|bin|bf|bf_trp]
    ├── vectors                     # (5) full pipeline in bf output mode (--mode bf_trp)
    │   ├── howde                   # with --mode bf_trp --split howde (input BFs for HowDeSBT cluster/build)
    │   └── sdsl                    # with --mode bf_trp --split sdsl
    ├── fof.txt                     # copy of input fof (--file)
    └── hash_window.vec             # info about partitioned bf (required for queries)
```

**CLI**
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
  --split STR                    split matrix in indidual files: [sdsl|howde] (only with --mode
                                 bf_trp) [default: none]
```

**Merge options**
* `--count-abundance-min INT`: An hard threshold, all k-mers with an abundance less than this parameter are discarded.
* `--merge-abundance-min INT/STR`: A soft threshold, all k-mers with an abundance between `count-abundance-min` and `merge-abundance-min` are considering rescue-able
* `--save-if INT`: If a k-mer is rescue-able, it is keep if it is solid (with an abundance greater than `merge-abundance-min`) in `save-if` other samples.
* `--recurrence-min INT`: All k-mers that do not occur in at least `recurrence-min` samples are discarded.

**Hash mode configuration (Bloom filters construction pipeline)**
* `--hasher STR`: hash function, for HowDe-SBT compatibility use `--hasher sabuhash`
* `--max-hash INT`: The requested Bloom filter size; kmtricks build partitioned Bloom filter (pBFs), each partition contains $\frac{max-hash}{nb-partitions} \text{(rounded up to 8)}$

**Output matrices**
* `--mode STR`:
  * k-mer mode
    * `ascii`: Each row corresponds to one k-mer and its counts, space-separated.
    * `bin`: Each row correponds to a k-mer and a count vector.
      * \<KType> [\<CType>, ..., \<CType>]
    * `pa`: Each row corresponds to a k-mer and a presence/absence bit-vector
      * \<KType> \<BitVector>
  * hash mode
    * `bf`: Each row corresponds to one presence/absence bit-vector corresponding to one hash value (given by the line index).
      * \<PaBitVector>
    * `bf_trp`: Each row represents
      * \<BfBitVector>

All matrix files have headers (in text format for `ascii` mode). When compression is enabled, headers remain uncompressed, compression is applied only on file content (even in `ascii` mode).

# A. Compute a minimizer repartition

**CLI example:**
```bash
kmtricks.py run --file fof.txt --run-dir ./km_minim_run --kmer-size 20 --nb-cores 8 --nb-partitions 50 --until repart
```
Repartition file is dumped at: `./km_minim_run/storage/partition_storage_gatb/minimRepart.minimRepart`

**Use the repartition:**
```cpp
#include <string>
#include <kmtricks/repartition.hpp>
#include <kmtricks/utilities.hpp>

using KType = km::selectK<32>::type;  // integral type that supports k up to 32
using CType = km::selectC<255>::type; // integral type that supports count up to 255

int main(int argc, char* argv[])
{
  const std::string path = "./km_minim_run/storage/partition_storage_gatb/minimRepart.minimRepart";
  km::RepartFile repartition(path);
  km::Kmer<KType> kmer("ACGTACGTACGTACGTACGT", true); // true = canonical
  km::Minimizer<KType> minimizer(&k, 10, true);       // true = discard low complexity minimizer
  auto partition = repartition.get(minimizer);
  return 0;
}
```

# B. Compute super-k-mers

**CLI example:**
```bash
kmtricks.py run --file fof.txt --run-dir ./superk_run --kmer-size 20 --nb-cores 8 --nb-partitions 50 --until superk
```

Super-k-mers are dumped at `./superk_run/storage/superk_partitions` with one directory per sample such as `${ID}.superk`. In these directories, each file corresponds to one partition: `superKparts.${part_id}`

**Read super-k-mers:**
```cpp
#include <string>
#include <kmtricks/skreader.hpp>
#include <kmtricks/utilities.hpp>

using KType = km::selectK<32>::type;  // integral type that supports k up to 32
using CType = km::selectC<255>::type; // integral type that supports count up to 255

int main(int argc, char* argv[])
{
  const std::string path = "./superk_run/storage/superk_partitions/D1.superk";
  const std::string pref = "superKparts.";
  const int nb_part      = 50;
  const size_t kmer_size = 20;
  
  km::SuperkStorage sk_store(path, pref, nb_part);
  
  km::SuperkReader<KType> sk_reader(&sk_store, kmer_size);
  km::Superk<KType> superkmer(kmer_size);
  
  int partition = 10; // read partition 10
  while (sk_reader.next_superk(partition, &superkmer))
    std::cout << superkmer.str_value() << std::endl;
  
  return 0;
}
```

# C. Count k-mers or hashes

## 1. k-mers mode

**CLI example:**
```bash
kmtricks.py run --file fof.txt --run-dir ./kmer_run --kmer-size 20 --nb-cores 8 --nb-partitions 50 --count-abundance-min 2 --until count --lz4
# --count-abundance-min is used only if abundances are not set in the fof
```

k-mers are dumped at `./kmer_run/storage/kmers_partitions` with one directory per partition. In these directories, each file corresponds to one sample: `${ID}.kmer` (with an optional `.lz4` extension if lz4 compression is enabled).

**Read a k-mer file:**

```cpp
#include <string>
#include <kmtricks/utilities.hpp>
#include <kmtricks/io.hpp>

using KType = km::selectK<32>::type;  // integral type that supports k up to 32
using CType = km::selectC<255>::type; // integral type that supports count up to 255

int main(int argc, char* argv[])
{
  // read k-mers from partition 1 from D1 sample
  const std::string path = "./km_run/storage/kmer_partitions/partition_0/D1.kmer.lz4";
  const std::string pref = "superKparts.";
  const int nb_part      = 50;
  const size_t kmer_size = 20;
  km::KmerFile<km::IN, KType, CType> kmer_file(path);
  CType count;

  km::Kmer<KType> kmer(false); // kmer is already canonical in KmerFile, check is not needed
  while (kmer_file.read(kmer, count))
    std::cout << kmer.str_value() << std::endl;
  
  return 0;
}
```

**Stream a count k-mer matrix:**

```cpp
#include <string>
#include <kmtricks/merger.hpp>
#include <kmtricks/utilities.hpp>
#include <kmtricks/io.hpp>

using KType = km::selectK<32>::type;  // integral type that supports k up to 32
using CType = km::selectC<255>::type; // integral type that supports count up to 255

int main(int argc, char* argv[])
{
  // stream partition 0
  const std::string path = "./kmer_run/storage/kmer_partitions/partition_0/partition0.fof";
  using F = km::KmerFile<km::IN, KType, CType>;
  const int merge_abundance_min = 5;
  const int recurrence_min = 1;
  const int kmer_size = 20;
  km::Merger<KType, CType, F> merger(fof, merge_abundance_min, recurrence_min, 0);
  
  // km::Kmer<KType> kmer(false); // optional
  while (!merger.end)
  {
    merger.next();
    if (merger.keep)
    {
      KType kmer_int = merger.m_khash; // kmer as integral type
      // kmer.set_kmer(merger.m_khash, kmer_size);
      for (auto& c: merger.counts)
        std::cout << std::to_string(c) << "\n";
    }
  }
  return 0;
}
```
## 2. hash counting mode

**CLI example:**
```bash
kmtricks.py run --file fof.txt --run-dir ./hash_run --kmer-size 20 --nb-cores 8 --nb-partitions 50 --count-abundance-min 2 --lz4 --hasher sabuhash --mode [bf | bf_trp] --hasher [sabuhash | xor] --max-hash 100000 --until count [--skip-merge]
```

With `--skip-merge`: hashes are represented by a bit-vector at `./kmer_run/storage/kmers_partitions/partition_${part_id}/${ID}.kmer.vec(.lz4)`  
Without `--skip-merge`: hashes and counts are dumped on disk at `./kmer_run/storage/kmers_partitions/partition_${part_id}/${ID}.kmer(.lz4)`

**Read a counted hashes file [without --skip-merge]:**

```cpp
#include <string>
#include <kmtricks/utilities.hpp>
#include <kmtricks/io.hpp>

using KType = km::selectK<32>::type;  // integral type that supports k up to 32
using CType = km::selectC<255>::type; // integral type that supports count up to 255

int main(int argc, char* argv[])
{
  const std::string path = "./kmer_run/storage/kmer_partitions/partition_0/D1.kmer.lz4";
  KmerFile<IN, uint64_t, CType> hash_file(path);

  uint64_t hash;
  CType count;

  while (hash_file.read(hash, count))
    std::cout << std::to_string(hash) << " " << std::to_string(count) << "\n";
}
```

**Read a bit vector file [with --skip-merge]:**

```cpp
#include <string>
#include <kmtricks/utilities.hpp>
#include <kmtricks/io.hpp>

using KType = km::selectK<32>::type;  // integral type that supports k up to 32
using CType = km::selectC<255>::type; // integral type that supports count up to 255

int main(int argc, char* argv[])
{
  const std::string path = "./kmer_run/storage/kmer_partitions/partition_0/D1.kmer.vec.lz4";
  BitVectorFile<IN> bitvector_file(path);
  vector<char> vec = bitvector_file.read();
  std::pair<uint64_t, uint64_t> window = bitvector_file.get_window();
}
```

# D. Build matrix

## 1. k-mer count matrix [k-mer mode]

**CLI example:**
```bash
kmtricks.py run --file fof.txt --run-dir ./count_matrix_run --kmer-size 20 --nb-cores 8 --nb-partitions 50 --count-abundance-min 2 --merge-abundance-min 4 --recurrence-min 1 --save-if 2 --mode [ascii | bin] --lz4 --until merge
```

One matrix is dumped per partition at `./count_matrix_run/storage/matrix/partition${part_id}/ascii_matrix${part_id}.mat` or at `./count_matrix_run/storage/matrix/partition${part_id}/count_matrix${part_id}.mat`.

**Read a count matrix:** 
```cpp
#include <string>
#include <vector>
#include <kmtricks/utilities.hpp>
#include <kmtricks/io.hpp>

using KType = km::selectK<32>::type;  // integral type that supports k up to 32
using CType = km::selectC<255>::type; // integral type that supports count up to 255

int main(int argc, char* argv[])
{
  km::CountMatrixFile<km::IN, KType, CType, km::matrix_t::ASCII> count_matrix (
    "./count_matrix_run/storage/matrix/partition_0/ascii_matrix0.mat"
  );

  // Or if it's a bin count matrix
  //km::CountMatrixFile<km::IN, KType, CType, km::matrix_t::BIN> count_matrix (
  //  "./count_matrix_run/storage/matrix/partition_0/count_matrix0.mat"
  //);
  
  km::Kmer<KType> kmer(false);
  std::vector<CType> counts(count_matrix.infos()->nb_counts);
  while (count_matrix.read(kmer, counts))
  {
    std::cout << kmer.str_value() << "\n";
    for (auto& c: counts)
      std::cout << " " << std::to_string(c);
    std::cout << "\n";
  }
}
```

## 2. k-mer presence/absence matrix [k-mer mode]

**CLI example:**
```bash
kmtricks.py run --file fof.txt --run-dir ./pa_matrix_run --kmer-size 20 --nb-cores 8 --nb-partitions 50 --count-abundance-min 2 --merge-abundance-min 4 --recurrence-min 1 --save-if 2 --mode pa --lz4 --until merge
```

One matrix is dumped per partition at `./pa_matrix_run/storage/matrix/partition${part_id}/pa_matrix${part_id}.mat`.

**Read a presence/absence matrix:**

```cpp
#include <string>
#include <vector>
#include <kmtricks/utilities.hpp>
#include <kmtricks/io.hpp>

using KType = km::selectK<32>::type;  // integral type that supports k up to 32
using CType = km::selectC<255>::type; // integral type that supports count up to 255

int main(int argc, char* argv[])
{
  km::PAMatrixFile<km::IN, KType> pa_matrix (
    "./pa_matrix_run/storage/matrix/partition_0/pa_matrix0.mat"
  );
  km::Kmer<KType> kmer(false);
  std::vector<char> pa_vectors(pa_matrix.infos()->size_in_bytes);
  while (pa.read(kmer, pa_vectors))
  {
    for (int file_id = 0; file_id<pa_matrix.infos()->bits_in_use; file_id++)
    {
      std::cout << kmer.str_value()
                << "in " << std::to_string(file_id)
                << std::boolalpha << BITCHECK(pa_vectors, file_id) << "\n"; 
    }
  }
}
```

## 3. Column-major Bloom filter matrix [hash counting mode]

**CLI example:**
```bash
kmtricks.py run --file fof.txt --run-dir ./cm_bf_matrix_run --kmer-size 20 --nb-cores 8 --nb-partitions 50 --count-abundance-min 2 --merge-abundance-min 4 --recurrence-min 1 --save-if 2 --mode bf --lz4 --until merge
```

One matrix is dumped per partition at `./cm_bf_matrix_run/storage/matrix/partition${part_id}/no_trp_bf${part_id}.mat`.

**Read a column-major bloom filter matrix:**

```cpp
#include <string>
#include <vector>
#include <kmtricks/utilities.hpp>
#include <kmtricks/io.hpp>

using KType = km::selectK<32>::type;  // integral type that supports k up to 32
using CType = km::selectC<255>::type; // integral type that supports count up to 255

int main(int argc, char* argv[])
{
  km::BitMatrixFile<km::IN, km::matrix_t::BF> bf_matrix (
    "./cm_bf_matrix_run/storage/matrix/partition_0/no_trp_bf0.mat"
  );
  std::pair<uint64_t, uint64_t> hash_window = bf_matrix.get_window();
  uint64_t lower_hash = hash_window.first;
  uint64_t upper_hash = hash_window.second;

  // here Bloom are column-major, so each line is a pres/abs vector corresponding to a hash value.
  vector<char> pa(bf_matrix.infos()->size_in_bytes)

  uint64_t current_hash = hash_window.first;
  while (bf_matrix.read(pa))
  {
    // at each loop 'pa' corresponds to a presence/absence vector for 'current_hash'
    current_hash++; // hashes are consecutive
  }
}
```

**Read a column-major bf matrix, tranpose it to obtain row-major bf (see next section) and write the transposed matrix into a kmtricks file:** 

```cpp
#include <string>
#include <vector>
#include <kmtricks/utilities.hpp>
#include <kmtricks/bitmatrix.hpp>
#include <kmtricks/io.hpp>

using KType = km::selectK<32>::type;  // integral type that supports k up to 32
using CType = km::selectC<255>::type; // integral type that supports count up to 255

int main(int argc, char* argv[])
{
  // km::matrix_t::BIT allows only to load/dump km::BitMatrix but format is the same as km::matrix_t::BF
  { // Read matrix, tranpose and write
    km::BitMatrixFile<km::IN, km::matrix_t::BIT> bf_matrix (
      "./cm_bf_matrix_run/storage/matrix/partition_0/no_trp_bf0.mat"
    );
    auto i = bf_matrix.infos()->nb_rows_use;
    auto j = bf_matrix.infos()->nb_cols_use/8;
    km::BitMatrix mat(i, j, true);
    bf_matrix.load(mat);

    km::BitMatrix* trp = mat.tranpose();
    auto i_trp = trp->get_nb_lines();
    auto j_trp = trp->get_nb_cols();

    int partition_id = 0;
    // add a compression layer (only for file content, headers are never compressed)
    int compress = 0; // without compression
    km::BitMatrixFile<km::OUT, km::matrix_t::BIT> out_matrix (
      "OutTranposeMatrix.mat", partition_id, i_trp, j_trp, compress
    );
    out_matrix.dump(*trp);
    delete trp;
  }
  { // Then reread use km::matrix_t::BF

  }
}
```

## 4. Row-major Bloom filter matrix (corresponding to a transposition of 3.) [hash counting mode]

**CLI example:**
```bash
kmtricks.py run --file fof.txt --run-dir ./rm_bf_matrix_run --kmer-size 20 --nb-cores 8 --nb-partitions 50 --count-abundance-min 2 --merge-abundance-min 4 --recurrence-min 1 --save-if 2 --mode bf_trp --lz4 --until merge
```

One matrix is dumped per partition at `./rm_bf_matrix_run/storage/matrix/partition${part_id}/trp_bf${part_id}.mat`.

**Read a row-major bf matrix**

```cpp
#include <string>
#include <vector>
#include <kmtricks/utilities.hpp>
#include <kmtricks/io.hpp>

using KType = km::selectK<32>::type;  // integral type that supports k up to 32
using CType = km::selectC<255>::type; // integral type that supports count up to 255

int main(int argc, char* argv[])
{
  km::BitMatrixFile<km::IN, km::matrix_t::BF> bf_matrix (
    "./rm_bf_matrix_run/storage/matrix/partition_0/trp_bf0.mat"
  );
  std::pair<uint64_t, uint64_t> hash_window = bf_matrix.get_window();
  uint64_t lower_hash = hash_window.first;
  uint64_t upper_hash = hash_window.second;

  std::vector<vector<char>> pBFs(bf_matrix.infos()->nb_rows,
                                 std::vector<char>(bf_matrix.infos()->nb_cols_use/8));
  int file_id = 0;
  
  while (bf_matrix.read(pBFs[file_id]))
  {
    file_id++;
  }
}
```

# E. Output [only in hash mode]

**CLI example:**
```bash
kmtricks.py run --file fof.txt --run-dir ./km_full_bf --kmer-size 20 --nb-cores 8 --nb-partitions 50 --count-abundance-min 2 --until count --lz4 --hasher sabuhash --mode bf_trp --hasher [sabuhash | xor] --max-hash 100000 --split [howde | sdsl]
```

* Parameters
  * `--split sdsl`: SDSL bit-vector compatibility
  * `--split howde`: HowDe-SBT compatibility

Individual Bloom filters are dumped at `./km_full_bf/storage/vectors/howde/${ID}.bf` or `./km_full_bf/storage/vectors/sdsl/${ID}.sdsl`.


# F. kmtricks and HowDe-SBT index

**Bloom filters construction with kmtricks**
```bash
kmtricks.py run --file fof.txt --run-dir ./km_index --kmer-size 20 --nb-cores 8 --nb-partitions 50 --count-abundance-min 2 --until count --lz4 --hasher sabuhash --mode bf_trp --hasher [sabuhash | xor] --max-hash 100000 --split [howde | sdsl]
```
**Index construction with HowDe-SBT**
```bash
cd ./km_index/storage/vectors/howde
ls *.bf > bf_list
km_howdesbt cluster bf_list
km_howdesbt build --howde bf_list
```

**Queries**
```bash
km_howdesbt queryKm --tree=bf_list.detbrief.sbt --repart=../../partition_storage_gatb/minimRepart.minimRepart --win=../../hash_window.vec query_file.fa > result.txt
```