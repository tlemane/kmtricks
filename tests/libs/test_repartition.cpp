#include <gtest/gtest.h>
#define _KM_LIB_INCLUDE_
#include <kmtricks/repartition.hpp>
#include <kmtricks/sequences.hpp>

using namespace km;
typedef uint64_t kt;

TEST(RepartFile, repartfile_build)
{
  string path = "./repartition_data/minimRepart.minimRepart";
  RepartFile repartition(path);
}

TEST(RepartFile, repartfile_get)
{
  string path = "./repartition_data/minimRepart.minimRepart";
  RepartFile repartition(path);

  Kmer<kt> kmer("CATACAGAGACAGCAGCAGA", true);
  Minimizer<kt> minimKmer(&kmer, 10, true);
  EXPECT_EQ(repartition.get(minimKmer.value()), 0);

  Superk<kt> superk("GAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  Minimizer<kt> minimSuperk(&superk, 10, true);
  EXPECT_EQ(repartition.get(minimSuperk.value()), 0);

  Kmer<kt> kmer1("TATTATATCTACTACCATCA", true);
  Minimizer<kt> minimKmer1(&kmer1, 10, true);
  EXPECT_EQ(repartition.get(minimKmer1.value()), 1);
}