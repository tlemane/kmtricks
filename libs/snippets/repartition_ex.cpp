#include "../kmtricks/repartition.hpp"
#include "../kmtricks/sequences.hpp"

using namespace std;
using namespace km;

int main(int argc, char* argv [])
{
  typedef uint64_t T;
  try
  {
    string path = argv[1];
    RepartFile repartition(path); // repartition file path, in kmtricks run: storage/partition_storage_gatb/minimRepart.minimRepart
    Kmer<T> k("GAGCAGCACAAACGAGACAC", true); // build a kmer
    Minimizer<T> m(&k, 10, true); // build a minimizer of size 10, with check validity = true
    cout << to_string(repartition.get(m.value())) << endl; // get partition
    cout << to_string(repartition.get(m)) << endl; // get partition
  }
  catch ( exception& e )
  {
    cout << e.what() << endl;
  }
  return 0;
}