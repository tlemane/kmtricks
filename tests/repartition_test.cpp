#include <string>
#include <gtest/gtest.h>
#include <kmtricks/kmer.hpp>
#include <kmtricks/repartition.hpp>

using namespace km;
TEST(repartition, repartition)
{
  Repartition repart("./data/repart_gatb/repartition.minimRepart", "");
  std::string k0 = "AATATACTATATAATATATATAGCGAGGGGG";
  std::string k1 = "AAAACGACGACCGCAACACGACGCCAGCAGA";
  std::string k2 = "AAGATATAATATATAAAATATATAGTGTCGT";
  std::string k3 = "AAAAAAAAAAAAAAAAAAAACGCGGCGAAAA";

  EXPECT_EQ(0, repart.get_partition(Kmer<32>(k0).minimizer(10).value()));
  EXPECT_EQ(1, repart.get_partition(Kmer<32>(k1).minimizer(10).value()));
  EXPECT_EQ(2, repart.get_partition(Kmer<32>(k2).minimizer(10).value()));
  EXPECT_EQ(3, repart.get_partition(Kmer<32>(k3).minimizer(10).value()));
}