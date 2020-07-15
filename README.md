# kmtricks

## kmtrics IOs

kmtricks is a **set of independent modules** designed for kmer counting given a set of raw read sets. 

**Input** is composed of a set of read sets in fasta or fastq format, gzipped or not.

**Final output** depends on the user usage. It may be

	* a matrix of kmer x abundance. $M_{i,j}$ is the abundance of kmer $i$ in the read set $j$
	* a matrix of kmer x presence or absence. $M_{i,j}$ is the presence (1) or absence (0) of kmer $i$ in the read set $j$
	* a matrix of hash_value x abundance. $M_{i,j}$ is the abundance of the hash_value $i$ in the read set $j$. 
	* a matrix of hash_value x presence or absence. $M_{i,j}$ is the presence (1) or absence (0) of the hash_value $i$ in the read set $j$. 

## kmtrics performances 

Compared to a usual pipeline as the one used by `HowDeSBT` using `JellyFish` for generating a bloom filter per input read set, `kmtricks` is 1.6 times faster. 

![perf_kmtricks_674](/Users/ppeterlo/workspace/kmtricks/doc/perf_kmtricks_674.png =250x)

Test realised with 674 RNA-seq experiments with 20 cores 100 GB RAM Intel(R) Xeon(R) 2.60GHz

List of IDs available [here](tests/kmtricks/experiment_list_674.txt).



