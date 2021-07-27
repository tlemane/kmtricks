#include <string>
#include <gtest/gtest.h>
#define private public
#include <gatb/gatb_core.hpp>
#include <kmtricks/utils.hpp>
#include <kmtricks/gatb/gatb_utils.hpp>


TEST(gatb_utils, copy_kmer32) // copy uint64_t
{
  typedef typename ::Kmer<32>::Type Type;

  std::string a = km::random_dna_seq(20);

  auto revc = [](char c) {return km::NToB[c];};

  km::Kmer<32> kmer; kmer.set_k(20);
  Type gkmer = Type::polynom(a.c_str(), 20, revc);

  km::copy_gatb_kmers(kmer, gkmer);
  EXPECT_EQ(kmer.to_string(), gkmer.toString(20));
}

TEST(gatb_utils, copy_kmer64) // copy __uint128_t
{
  typedef typename ::Kmer<64>::Type Type2;

  std::string b = km::random_dna_seq(40);

  auto revc = [](char c) {return km::NToB[c];};

  km::Kmer<64> kmer2; kmer2.set_k(40);
  Type2 gkmer2 = Type2::polynom(b.c_str(), 40, revc);

  km::copy_gatb_kmers(kmer2, gkmer2);
  EXPECT_EQ(kmer2.to_string(), gkmer2.toString(40));
}

TEST(gatb_utils, copy_kmerX) // copy uint64_t[]
{
  typedef typename ::Kmer<92>::Type Type3;

  std::string c = km::random_dna_seq(90);

  auto revc = [](char c) {return km::NToB[c];};

  km::Kmer<92> kmer3; kmer3.set_k(90);
  Type3 gkmer3 = Type3::polynom(c.c_str(), 90, revc);

  km::copy_gatb_kmers(kmer3, gkmer3);
  EXPECT_EQ(kmer3.to_string(), gkmer3.toString(90));
}
