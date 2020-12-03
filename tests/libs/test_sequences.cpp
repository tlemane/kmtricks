#include <gtest/gtest.h>
#define _KM_LIB_INCLUDE_
#include <kmtricks/sequences.hpp>

using namespace km;
typedef uint64_t kt;

TEST(kmer, kmer_build_from_string)
{
  Kmer<kt> kmer("ACGTACGT", false);
  EXPECT_EQ(kmer.value(), 0x1E1E);
  EXPECT_EQ(kmer.str_value(), "ACGTACGT");
}

TEST(kmer, kmer_build_from_value)
{
  Kmer<kt> kmer(0x1E1E, 8, false);
  EXPECT_EQ(kmer.value(), 0x1E1E);
  EXPECT_EQ(kmer.str_value(), "ACGTACGT");
} 

TEST(kmer, kmer_rev_comp)
{
  Kmer<kt> kmer("ACGTACTT", false);
  EXPECT_EQ(kmer.str_rev_comp(), "AAGTACGT");
  EXPECT_EQ(kmer.rev_comp(), 0x0E1E);
}

TEST(kmer, kmer_hash)
{
  Kmer<kt> kmer("ACGTACGT", false);
  EXPECT_EQ(kmer.hash(), 0x4BC4D2729806CDF8);
}

template<typename K>
class TestCustomHasher : public Hasher<K>
{
public:
  TestCustomHasher() = default;
  ~TestCustomHasher() = default;
  uint64_t operator()(K data, uint64_t seed)
  {
    return 0x1234;
  }
};

TEST(kmer, kmer_custom_hash)
{
  TestCustomHasher<kt> *hasher = new TestCustomHasher<kt>();
  Kmer<kt> kmer("ACGTACGT", false);
  EXPECT_EQ(kmer.hash(), 0x4BC4D2729806CDF8);
  EXPECT_EQ(kmer.hash(hasher), 0x1234);
  EXPECT_EQ(kmer.hash(), 0x4BC4D2729806CDF8);
  kmer.set_hasher(hasher);
  EXPECT_EQ(kmer.hash(), 0x1234);
  kmer.set_default_hasher();
  EXPECT_EQ(kmer.hash(), 0x4BC4D2729806CDF8);
}

TEST(kmer, kmer_custom_encoding)
{
  uchar map[4] = {'T','A','C','G'};
  Code<kt> my_code(map);
  Kmer<kt> kmer("ACGTACGT", false, &my_code); // 0110110001101100
  EXPECT_EQ(kmer.value(), 0x6C6C);
}

TEST(kmer, kmer_operator)
{
  Kmer<kt> kmer("ACGTACGT", false);
  EXPECT_LT(kmer, 0xFFFF);
  EXPECT_GT(kmer, 0xFF);
  EXPECT_NE(kmer, 0xFF);
  EXPECT_EQ(kmer, 0x1E1E);

  EXPECT_LT(kmer, "AGCGTACG");
  EXPECT_GT(kmer, "AAGTACGT");
  EXPECT_NE(kmer, "ACGTACGG");
  EXPECT_EQ(kmer, "ACGTACGT");
}

TEST(kmer, kmer_canonical)
{
  Kmer<kt> kmer("ACGTACTT", true);
  EXPECT_EQ(kmer.str_value(), "AAGTACGT");
  EXPECT_EQ(kmer.value(), 0x0E1E);

  Kmer<kt> kmer2("ACGTACTT", false);
  kmer2.use_canonical();
  EXPECT_EQ(kmer2.str_value(), "AAGTACGT");
  EXPECT_EQ(kmer2.value(), 0x0E1E);
  EXPECT_EQ(kmer2.rev_comp(), 0x1E1A);
  EXPECT_EQ(kmer2.str_rev_comp(), "ACGTACTT");
}

TEST(kmer, kmer_set)
{
  Kmer<kt> kmer("ACGTACGT", false);
  EXPECT_EQ(kmer.str_value(), "ACGTACGT");
  kmer.set_kmer("ACCCTTTA");
  EXPECT_EQ(kmer.value(), 0x15A8);
  EXPECT_EQ(kmer.str_value(), "ACCCTTTA");
}

TEST(superk, superk_build_from_string)
{
  Superk<kt> superk("GAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  EXPECT_EQ(superk.str_value(), "GAGCAGCACAAACGAGACACAAAAAAAGAG");
  EXPECT_EQ(superk.value()[0], 0xCD);
  EXPECT_EQ(superk.value()[1], 0x34);
  EXPECT_EQ(superk.value()[2], 0x40);
  EXPECT_EQ(superk.value()[3], 0x73);
  EXPECT_EQ(superk.value()[4], 0x11);
  EXPECT_EQ(superk.value()[5], 0x00);
  EXPECT_EQ(superk.value()[6], 0x03);
  EXPECT_EQ(superk.value()[7], 0x30);
}

TEST(superk, superk_build_from_buffer)
{
  uchar buf[8] = {0xCD, 0x34, 0x40, 0x73, 0x11, 0x00, 0x03, 0x30};
  Superk<kt> superk(buf, 30, 20);
  EXPECT_EQ(superk.str_value(), "GAGCAGCACAAACGAGACACAAAAAAAGAG");
}

TEST(superk, superk_build_from_gatb_format)
{
  uchar buf[8] = {0x11, 0x73, 0x40, 0x34, 0xCD, 0x00, 0xC0, 0x0C};
  Superk<kt> superk(buf, 30, 20, true);
  EXPECT_EQ(superk.str_value(), "GAGCAGCACAAACGAGACACAAAAAAAGAG");
}

TEST(superk, superk_set_from_string)
{
  Superk<kt> superk("GAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  EXPECT_EQ(superk.str_value(), "GAGCAGCACAAACGAGACACAAAAAAAGAG");
  superk.set_superk("CATACAGAGACAGCAGCAGAGCA");
  EXPECT_EQ(superk.str_value(), "CATACAGAGACAGCAGCAGAGCA");
}

TEST(superk, superk_set_from_buffer)
{
  Superk<kt> superk("GAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  EXPECT_EQ(superk.str_value(), "GAGCAGCACAAACGAGACACAAAAAAAGAG");
  uchar buf[6] = {0x48, 0x4C , 0xC4, 0xD3, 0x4C, 0xD0};
  superk.set_superk(buf, 23, 20);
  EXPECT_EQ(superk.str_value(), "CATACAGAGACAGCAGCAGAGCA");
}

TEST(superk, superk_set_from_gatb_format)
{
  Superk<kt> superk("CATACAGAGACAGCAGCAGAGCA", 20);
  EXPECT_EQ(superk.str_value(), "CATACAGAGACAGCAGCAGAGCA");
  uchar buf[8] = {0x11, 0x73, 0x40, 0x34, 0xCD, 0x00, 0xC0, 0x0C};
  superk.set_superk(buf, 30, 20, true);
  EXPECT_EQ(superk.str_value(), "GAGCAGCACAAACGAGACACAAAAAAAGAG");
}

TEST(superk, superk_get_kmer)
{
  Superk<kt> superk("GAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  EXPECT_EQ(superk.get_kmer(0, false).value(), 0xCD34407311);
  EXPECT_EQ(superk.get_kmer(1, false).value(), 0x34D101CC44);
}

TEST(superk, superk_nb_kmers)
{
  Superk<kt> superk("GAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  EXPECT_EQ(superk.nb_kmers(), 11);
}

TEST(superk, superk_operator)
{
  Superk<kt> superk("GAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  Superk<kt> superk2("GAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  Superk<kt> superk3("AAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  Superk<kt> superk4("TAGCAGCACAAACGAGACACAAAAAAAGAG", 20);

  EXPECT_EQ(superk, superk2);
  EXPECT_LT(superk, superk4);
  EXPECT_GT(superk, superk3);
  EXPECT_NE(superk, superk3);

  EXPECT_EQ(superk, "GAGCAGCACAAACGAGACACAAAAAAAGAG");
  EXPECT_LT(superk, "TAGCAGCACAAACGAGACACAAAAAAAGAG");
  EXPECT_GT(superk, "AAGCAGCACAAACGAGACACAAAAAAAGAG");
  EXPECT_NE(superk, "AAGCAGCACAAACGAGACACAAAAAAAGAG");
}

TEST(minimizer, minimizer_build_from_kmer)
{
  Kmer<kt> kmer("GAGCAGCACAAACGAGACAC", true);
  Minimizer<kt> minim(&kmer, 10, false);
  EXPECT_EQ(minim.value(), 0x1CC4);
  EXPECT_EQ(minim.str_value(), "AAACGAGACA"); // not valid minimizer
  Minimizer<kt> minim2(&kmer, 10, true); // check validity
  EXPECT_EQ(minim2.value(), 0x7311);
  EXPECT_EQ(minim2.str_value(), "AACGAGACAC");
}

TEST(minimizer, minimizer_build_from_superk)
{
  Superk<kt> superk("GAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  Minimizer<kt> minim(&superk, 10, false);
  EXPECT_EQ(minim.value(), 0x1CC4);
  EXPECT_EQ(minim.str_value(), "AAACGAGACA");
}

TEST(minimizer, minimizer_set_from_kmer)
{
  Kmer<kt> kmer("GAGCAGCACAAACGAGACAC", true);
  Minimizer<kt> minim(10);
  minim.set_kmer(&kmer, 10, true);
  EXPECT_EQ(minim.value(), 0x7311);
}

TEST(minimizer, minimizer_set_from_superk)
{
  Superk<kt> superk("GAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  Minimizer<kt> minim(10);
  minim.set_superk(&superk, 10, true);
  EXPECT_EQ(minim.value(), 0x7311);
}

TEST(minimizer, minimizer_operator)
{
  Kmer<kt> kmer("GAGCAGCACAAACGAGACAC", true);
  Minimizer<kt> minim(&kmer, 10, false);

  EXPECT_LT(minim, 0xFFFF);
  EXPECT_GT(minim, 0xFF);
  EXPECT_NE(minim, 0xFF);
  EXPECT_EQ(minim, 0x1CC4);

  EXPECT_LT(minim, "AGCGTACGAA");
  EXPECT_GT(minim, "AAAAACGTAA");
  EXPECT_NE(minim, "ACGTACGGAA");
  EXPECT_EQ(minim, "AAACGAGACA");
}

TEST(minimizer, minimizer_default)
{
  Kmer<kt> kmer("AAAAAAAAAAAAAAAAAAAAAAA", true);
  Minimizer<kt> minim(&kmer, 10, true);
  EXPECT_EQ(minim.value(), DEFAULT_MINIMIZER);

  minim.set_default(0x1234);
  EXPECT_EQ(minim.value(), 0x1234);

  minim.set_default("ACGTACGTAA");
  EXPECT_EQ(minim.value(), 0x01E1E0);

  minim.set_default();
  EXPECT_EQ(minim.value(), DEFAULT_MINIMIZER);
}