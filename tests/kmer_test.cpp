#include <string>
#include <gtest/gtest.h>
#define private public
#include <kmtricks/kmer.hpp>
#include <kmtricks/utils.hpp>

using namespace km;

TEST(kmer, set_from_str)
{
  std::string a = km::random_dna_seq(20);
  std::string b = km::random_dna_seq(40);
  std::string c = km::random_dna_seq(90);
  std::string d = km::random_dna_seq(32);

  km::Kmer<32> kmer(a);
  EXPECT_EQ(kmer.name(), "Kmer<32> - uint64_t");
  EXPECT_EQ(a , kmer.to_string());

  km::Kmer<64> kmer2(b);
  EXPECT_EQ(kmer2.name(), "Kmer<64> - __uint128_t");
  EXPECT_EQ(b , kmer2.to_string());

  km::Kmer<92> kmer3(c);
  EXPECT_EQ(kmer3.name(), "Kmer<92> - uint64_t[3]");
  EXPECT_EQ(c , kmer3.to_string());

  km::Kmer<32> kmer4(d);
  EXPECT_EQ(d , kmer4.to_string());
}

TEST(kmer, base)
{
  std::string a = km::random_dna_seq(20);
  std::string b = km::random_dna_seq(40);
  std::string c = km::random_dna_seq(90);

  km::Kmer<32> kmer(a);
  for (size_t i=0; i<a.size(); i++)
    EXPECT_EQ(a[i], kmer.at(i));

  km::Kmer<64> kmer2(b);
  for (size_t i=0; i<b.size(); i++)
    EXPECT_EQ(b[i], kmer2.at(i));

  km::Kmer<92> kmer3(c);
  for (size_t i=0; i<c.size(); i++)
    EXPECT_EQ(c[i], kmer3.at(i));
}

TEST(kmer, rev_comp)
{
  std::string a = km::random_dna_seq(20); std::string a2 = km::str_rev_comp(a);
  std::string b = km::random_dna_seq(40); std::string b2 = km::str_rev_comp(b);
  std::string c = km::random_dna_seq(90); std::string c2 = km::str_rev_comp(c);

  km::Kmer<32> kmer(a);
  km::Kmer<32> rkmer = kmer.rev_comp();
  EXPECT_EQ(a2, rkmer.to_string());

  km::Kmer<64> kmer2(b);
  km::Kmer<64> rkmer2 = kmer2.rev_comp();
  EXPECT_EQ(b2, rkmer2.to_string());

  km::Kmer<92> kmer3(c);
  km::Kmer<92> rkmer3 = kmer3.rev_comp();
  EXPECT_EQ(c2, rkmer3.to_string());
}

TEST(kmer, canonical)
{
  std::string a = "AAAAAAACCCCCCC";
  std::string b = "CGCCCCCCCCCCCT";
  std::string c = "AGGGGGGGGGGGCG";

  km::Kmer<32> kmer(a);
  EXPECT_EQ(kmer.canonical().to_string(), a);

  km::Kmer<32> kmer2(b);
  EXPECT_EQ(kmer2.canonical().to_string(), c);
  EXPECT_NE(kmer2.canonical().to_string(), b);
}

TEST(kmer, operator)
{
  std::string a = "AAAAAAACCCCCCC";
  std::string b = "AAAAAAACCCCCCT";

  std::string c = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCC";
  std::string d = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCT";

  std::string e = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCC";
  std::string f = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCT";

  km::Kmer<32> kmera(a); km::Kmer<32> kmerb(b);
  EXPECT_TRUE(kmera < kmerb);
  EXPECT_FALSE(kmera > kmerb);
  EXPECT_FALSE(kmera == kmerb);
  EXPECT_TRUE(kmera != kmerb);
  EXPECT_TRUE(kmera == kmera);

  km::Kmer<64> kmerc(c); km::Kmer<64> kmerd(d);
  EXPECT_TRUE(kmerc < kmerd);
  EXPECT_FALSE(kmerc > kmerd);
  EXPECT_FALSE(kmerc == kmerd);
  EXPECT_TRUE(kmerc != kmerd);
  EXPECT_TRUE(kmerc == kmerc);

  km::Kmer<96> kmere(e); km::Kmer<96> kmerf(f);
  EXPECT_TRUE(kmere < kmerf);
  EXPECT_FALSE(kmere > kmerf);
  EXPECT_FALSE(kmere == kmerf);
  EXPECT_TRUE(kmere != kmerf);
  EXPECT_TRUE(kmere == kmere);
}

TEST(kmer, minimizer)
{
  {
    std::string a = "ACGAGCAATACGA";
    Kmer<32> kmer(a);
    auto v = kmer.mmers(4);
    EXPECT_EQ("ACGA", v[0].to_string());
    EXPECT_EQ("CGAG", v[1].to_string());
    EXPECT_EQ("GAGC", v[2].to_string());
    EXPECT_EQ("AGCA", v[3].to_string());
    EXPECT_EQ("GCAA", v[4].to_string());
    EXPECT_EQ("CAAT", v[5].to_string());
    EXPECT_EQ("AATA", v[6].to_string());
    EXPECT_EQ("ATAC", v[7].to_string());
    EXPECT_EQ("TACG", v[8].to_string());
    EXPECT_EQ("ACGA", v[9].to_string());

    Mmer m = kmer.minimizer(4);
    EXPECT_EQ(m.to_string(), "AATA");
  }
  {
    std::string a = "ACGAGCAATACGA";
    Kmer<32> kmer(a);
    auto v = kmer.mmers(4);
    EXPECT_EQ("ACGA", v[0].to_string());
    EXPECT_EQ("CGAG", v[1].to_string());
    EXPECT_EQ("GAGC", v[2].to_string());
    EXPECT_EQ("AGCA", v[3].to_string());
    EXPECT_EQ("GCAA", v[4].to_string());
    EXPECT_EQ("CAAT", v[5].to_string());
    EXPECT_EQ("AATA", v[6].to_string());
    EXPECT_EQ("ATAC", v[7].to_string());
    EXPECT_EQ("TACG", v[8].to_string());
    EXPECT_EQ("ACGA", v[9].to_string());

    Mmer m = kmer.minimizer(4);
    EXPECT_EQ(m.to_string(), "AATA");
  }
}