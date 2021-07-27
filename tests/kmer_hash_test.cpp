#include <gtest/gtest.h>
#define WITH_XXHASH
#include <kmtricks/kmer_hash.hpp>
#include <kmtricks/utils.hpp>

using namespace km;

TEST(kmer, kmer_hash_template)
{
  EXPECT_EQ(KmerHashers<0>::name(), "KmerHashers<0> - Folly hash");
  // uint64_t spe
  EXPECT_EQ(KmerHashers<0>::Hasher<32>::name(), "KmerHashers<0>::Hasher<32>");
  EXPECT_EQ(KmerHashers<0>::WinHasher<32>::name(), "KmerHashers<0>::WinHasher<32>");

  // __uint128_t spe
  EXPECT_EQ(KmerHashers<0>::Hasher<64>::name(), "KmerHashers<0>::Hasher<64>");
  EXPECT_EQ(KmerHashers<0>::WinHasher<64>::name(), "KmerHashers<0>::WinHasher<64>");

  // Generic for K > 64
  EXPECT_EQ(KmerHashers<0>::Hasher<96>::name(), "KmerHashers<0>::Hasher<MAX_K=96>");
  EXPECT_EQ(KmerHashers<0>::Hasher<128>::name(), "KmerHashers<0>::Hasher<MAX_K=128>");
  EXPECT_EQ(KmerHashers<0>::WinHasher<96>::name(), "KmerHashers<0>::WinHasher<MAX_K=96>");
  EXPECT_EQ(KmerHashers<0>::WinHasher<128>::name(), "KmerHashers<0>::WinHasher<MAX_K=128>");

  EXPECT_EQ(KmerHashers<1>::name(), "KmerHashers<1> - XXHASH");
  // uint64_t spe
  EXPECT_EQ(KmerHashers<1>::Hasher<32>::name(), "KmerHashers<1>::Hasher<32>");
  EXPECT_EQ(KmerHashers<1>::WinHasher<32>::name(), "KmerHashers<1>::WinHasher<32>");

  // __uint128_t spe
  EXPECT_EQ(KmerHashers<1>::Hasher<64>::name(), "KmerHashers<1>::Hasher<64>");
  EXPECT_EQ(KmerHashers<1>::WinHasher<64>::name(), "KmerHashers<1>::WinHasher<64>");

  // Generic for K > 64
  EXPECT_EQ(KmerHashers<1>::Hasher<96>::name(), "KmerHashers<1>::Hasher<MAX_K=96>");
  EXPECT_EQ(KmerHashers<1>::Hasher<128>::name(), "KmerHashers<1>::Hasher<MAX_K=128>");
  EXPECT_EQ(KmerHashers<1>::WinHasher<96>::name(), "KmerHashers<1>::WinHasher<MAX_K=96>");
  EXPECT_EQ(KmerHashers<1>::WinHasher<128>::name(), "KmerHashers<1>::WinHasher<MAX_K=128>");
}

TEST(kmer, folly_hash)
{
  using Folly = KmerHashers<0>; // use folly hash
  {
    std::string k1 = random_dna_seq(20);
    Kmer<32> ka(k1);
    Kmer<32> kb(k1);
    Folly::Hasher<32> hasher;
    Folly::WinHasher<32> winhasher(0, 1000);
    EXPECT_EQ(hasher(ka), hasher(kb));
    EXPECT_EQ(winhasher(ka), winhasher(kb));
    EXPECT_LT(winhasher(ka), 1001);
  }
  {
    std::string k1 = random_dna_seq(40);
    Kmer<64> ka(k1);
    Kmer<64> kb(k1);
    Folly::Hasher<64> hasher;
    Folly::WinHasher<64> winhasher(0, 1000);
    EXPECT_EQ(hasher(ka), hasher(kb));
    EXPECT_EQ(winhasher(ka), winhasher(kb));
    EXPECT_LT(winhasher(ka), 1001);
  }
  {
    std::string k1 = random_dna_seq(90);
    Kmer<96> ka(k1);
    Kmer<96> kb(k1);
    Folly::Hasher<96> hasher;
    Folly::WinHasher<96> winhasher(0, 1000);
    EXPECT_EQ(hasher(ka), hasher(kb));
    EXPECT_EQ(winhasher(ka), winhasher(kb));
    EXPECT_LT(winhasher(ka), 1001);
  }
}

TEST(kmer, xxhash_hash)
{
  using XXHash = KmerHashers<1>; // use xxhash
  {
    std::string k1 = random_dna_seq(20);
    Kmer<32> ka(k1);
    Kmer<32> kb(k1);
    XXHash::Hasher<32> hasher;
    XXHash::WinHasher<32> winhasher(0, 1000);
    EXPECT_EQ(hasher(ka), hasher(kb));
    EXPECT_EQ(winhasher(ka), winhasher(kb));
    EXPECT_LT(winhasher(ka), 1001);
  }
  {
    std::string k1 = random_dna_seq(20);
    Kmer<64> ka(k1);
    Kmer<64> kb(k1);
    XXHash::Hasher<64> hasher;
    XXHash::WinHasher<64> winhasher(0, 1000);
    EXPECT_EQ(hasher(ka), hasher(kb));
    EXPECT_EQ(winhasher(ka), winhasher(kb));
    EXPECT_LT(winhasher(ka), 1001);
  }
  {
    std::string k1 = random_dna_seq(20);
    Kmer<96> ka(k1);
    Kmer<96> kb(k1);
    XXHash::Hasher<96> hasher;
    XXHash::WinHasher<96> winhasher(0, 1000);
    EXPECT_EQ(hasher(ka), hasher(kb));
    EXPECT_EQ(winhasher(ka), winhasher(kb));
    EXPECT_LT(winhasher(ka), 1001);
  }
}