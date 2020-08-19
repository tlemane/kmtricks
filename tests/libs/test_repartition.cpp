#include <criterion/criterion.h>
#include "repartition.hpp"
#include "sequences.hpp"

using namespace km;
typedef uint64_t kt;

Test(RepartFile, repartfile_build)
{
  string path = "./repartition_data/minimRepart.minimRepart";
  RepartFile repartition(path);
}

Test(RepartFile, repartfile_get)
{
  string path = "./repartition_data/minimRepart.minimRepart";
  RepartFile repartition(path);

  Kmer<kt> kmer("CATACAGAGACAGCAGCAGA", true);
  Minimizer<kt> minimKmer(&kmer, 10, true);
  cr_assert(repartition.get(minimKmer.value()) == 0);

  Superk<kt> superk("GAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  Minimizer<kt> minimSuperk(&superk, 10, true);
  cr_assert(repartition.get(minimSuperk.value()) == 0);

  Kmer<kt> kmer1("TATTATATCTACTACCATCA", true);
  Minimizer<kt> minimKmer1(&kmer1, 10, true);
  cr_assert(repartition.get(minimKmer1.value()) == 1);
}