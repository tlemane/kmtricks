#include <gtest/gtest.h>
#include <kmtricks/io/kmer_file.hpp>
#include <kmtricks/utils.hpp>

using namespace km;

TEST(kmer_file, KmerWriter)
{
  KmerWriter kw("tests_tmp/k1.kmer", 21, 1, 1, 2, false);
  EXPECT_EQ(kw.infos().kmer_size, 21);
  EXPECT_EQ(kw.infos().kmer_slots, (21+31)/32);
  EXPECT_EQ(kw.infos().count_slots, 1);
  EXPECT_EQ(kw.infos().id, 1);
  EXPECT_EQ(kw.infos().partition, 2);
  EXPECT_EQ(kw.infos().compressed, false);
}

TEST(kmer_file, KmerReader)
{
  KmerReader kr("tests_tmp/k1.kmer");
  EXPECT_EQ(kr.infos().kmer_size, 21);
  EXPECT_EQ(kr.infos().kmer_slots, (21+31)/32);
  EXPECT_EQ(kr.infos().count_slots, 1);
  EXPECT_EQ(kr.infos().id, 1);
  EXPECT_EQ(kr.infos().partition, 2);
  EXPECT_EQ(kr.infos().compressed, false);
}

TEST(kmer_file, KmerWriteRead)
{
  std::vector<std::string> str_kmers(10000);
  {
    KmerWriter kw("tests_tmp/k2.kmer", 21, 1, 1, 2, false);
    KmerWriter kw2("tests_tmp/k2.kmer.lz4", 21, 1, 1, 2, true);
    for (size_t i=0; i<str_kmers.size(); i++)
    {
      str_kmers[i] = random_dna_seq(21);
      Kmer<32> kmer(str_kmers[i]);
      kw.write<32, 255>(kmer, 42);
      kw2.write<32, 255>(kmer, 42);
    }
  }
  {
    KmerReader kr("tests_tmp/k2.kmer");
    KmerReader kr2("tests_tmp/k2.kmer.lz4");
    Kmer<32> kmer; kmer.set_k(kr.infos().kmer_size);
    uint8_t c = 0;
    for (size_t i=0; i<str_kmers.size(); i++)
    {
      kr.read<32, 255>(kmer, c);
      EXPECT_EQ(kmer.to_string(), str_kmers[i]);
      EXPECT_EQ(c, 42);
      kr2.read<32, 255>(kmer, c);
      EXPECT_EQ(kmer.to_string(), str_kmers[i]);
      EXPECT_EQ(c, 42);
    }
  }
  {
    std::sort(str_kmers.begin(), str_kmers.end());
    KmerWriter kw("tests_tmp/k3.kmer.lz4", 21, 1, 1, 2, true);
    for (size_t i=0; i<str_kmers.size(); i++)
    {
      Kmer<32> kmer(str_kmers[i]);
      kw.write<32, 255>(kmer, 42);
    }
  }
  {
    std::ofstream out("tests_tmp/k3.kmer.csv");
    KmerReader("tests_tmp/k3.kmer.lz4").write_as_text<32, 255>(out);
  }
}
