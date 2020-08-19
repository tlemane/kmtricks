#include "kmtricks/skreader.hpp"

using namespace km;

int main(int argc, char* argv[])
{
  // Super-k-mers storage: directory, super-k-mers files prefix
  string path = argv[1];
  string prefix = "superKparts.";
  int nbpart = 4;
  SuperkStorage store(path, prefix, nbpart);

  // Super-k-mers reader: storage, k-mer size
  size_t kmer_size = 20;
  SuperkReader<uint64_t> reader(&store, kmer_size);

  // superk to insert each super-k-mers read
  Superk<uint64_t> superk(kmer_size);

  int partition = 0;
  while ( reader.next_superk(partition, &superk))
  {
    cout << to_string(superk.size()) << " " << superk.str_value() << endl;
  }
}