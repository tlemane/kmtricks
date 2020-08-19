#include <criterion/criterion.h>
#include "sequences.hpp"

using namespace km;
typedef uint64_t kt;

Test(kmer, kmer_build_from_string)
{
  Kmer<kt> kmer("ACGTACGT", false);
  cr_assert(kmer.value() == 0x1E1E);
  cr_assert(kmer.str_value() == "ACGTACGT");
}

Test(kmer, kmer_build_from_value)
{
  Kmer<kt> kmer(0x1E1E, 8, false);
  cr_assert(kmer.value() == 0x1E1E);
  cr_assert(kmer.str_value() == "ACGTACGT");
}

Test(kmer, kmer_rev_comp)
{
  Kmer<kt> kmer("ACGTACTT", false);
  cr_assert(kmer.str_rev_comp() == "AAGTACGT");
  cr_assert(kmer.rev_comp() == 0x0E1E);
}

Test(kmer, kmer_hash)
{
  Kmer<kt> kmer("ACGTACGT", false);
  cr_assert(kmer.hash() == 0x4BC4D2729806CDF8);
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

Test(kmer, kmer_custom_hash)
{
  TestCustomHasher<kt> *hasher = new TestCustomHasher<kt>();
  Kmer<kt> kmer("ACGTACGT", false);
  cr_assert(kmer.hash() == 0x4BC4D2729806CDF8);
  cr_assert(kmer.hash(hasher) == 0x1234);
  cr_assert(kmer.hash() == 0x4BC4D2729806CDF8);
  kmer.set_hasher(hasher);
  cr_assert(kmer.hash() == 0x1234);
  kmer.set_default_hasher();
  cr_assert(kmer.hash() == 0x4BC4D2729806CDF8);
}

Test(kmer, kmer_custom_encoding)
{
  uchar map[4] = {'T','A','C','G'};
  Code<kt> my_code(map);
  Kmer<kt> kmer("ACGTACGT", false, &my_code); // 0110110001101100
  cr_assert(kmer.value() == 0x6C6C);
}

Test(kmer, kmer_operator)
{
  Kmer<kt> kmer("ACGTACGT", false);
  cr_assert(kmer < 0xFFFF);
  cr_assert(kmer > 0xFF);
  cr_assert(kmer != 0xFF);
  cr_assert(kmer == 0x1E1E);

  cr_assert(kmer < "AGCGTACG");
  cr_assert(kmer > "AAGTACGT");
  cr_assert(kmer != "ACGTACGG");
  cr_assert(kmer == "ACGTACGT");
}

Test(kmer, kmer_canonical)
{
  Kmer<kt> kmer("ACGTACTT", true);
  cr_assert(kmer.str_value() == "AAGTACGT");
  cr_assert(kmer.value() == 0x0E1E);

  Kmer<kt> kmer2("ACGTACTT", false);
  kmer2.use_canonical();
  cr_assert(kmer2.str_value() == "AAGTACGT");
  cr_assert(kmer2.value() == 0x0E1E);
  cr_assert(kmer2.rev_comp() == 0x1E1A);
  cr_assert(kmer2.str_rev_comp() == "ACGTACTT");
}

Test(kmer, kmer_set)
{
  Kmer<kt> kmer("ACGTACGT", false);
  cr_assert(kmer.str_value() ==  "ACGTACGT");
  kmer.set_kmer("ACCCTTTA");
  cr_assert(kmer.value() == 0x15A8);
  cr_assert(kmer.str_value() == "ACCCTTTA");
}

Test(superk, superk_build_from_string)
{
  Superk<kt> superk("GAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  cr_assert(superk.str_value() == "GAGCAGCACAAACGAGACACAAAAAAAGAG");
  cr_assert(superk.value()[0] == 0xCD);
  cr_assert(superk.value()[1] == 0x34);
  cr_assert(superk.value()[2] == 0x40);
  cr_assert(superk.value()[3] == 0x73);
  cr_assert(superk.value()[4] == 0x11);
  cr_assert(superk.value()[5] == 0x00);
  cr_assert(superk.value()[6] == 0x03);
  cr_assert(superk.value()[7] == 0x30);
}

Test(superk, superk_build_from_buffer)
{
  uchar buf[8] = {0xCD, 0x34, 0x40, 0x73, 0x11, 0x00, 0x03, 0x30};
  Superk<kt> superk(buf, 30, 20);
  cr_assert(superk.str_value() == "GAGCAGCACAAACGAGACACAAAAAAAGAG");
}

Test(superk, superk_build_from_gatb_format)
{
  uchar buf[8] = {0x11, 0x73, 0x40, 0x34, 0xCD, 0x00, 0xC0, 0x0C};
  Superk<kt> superk(buf, 30, 20, true);
  cr_assert(superk.str_value() == "GAGCAGCACAAACGAGACACAAAAAAAGAG");
}

Test(superk, superk_set_from_string)
{
  Superk<kt> superk("GAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  cr_assert(superk.str_value() == "GAGCAGCACAAACGAGACACAAAAAAAGAG");
  superk.set_superk("CATACAGAGACAGCAGCAGAGCA");
  cr_assert(superk.str_value() == "CATACAGAGACAGCAGCAGAGCA");
}

Test(superk, superk_set_from_buffer)
{
  Superk<kt> superk("GAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  cr_assert(superk.str_value() == "GAGCAGCACAAACGAGACACAAAAAAAGAG");
  uchar buf[6] = {0x48, 0x4C , 0xC4, 0xD3, 0x4C, 0xD0};
  superk.set_superk(buf, 23, 20);
  cr_assert(superk.str_value() == "CATACAGAGACAGCAGCAGAGCA");
}

Test(superk, superk_set_from_gatb_format)
{
  Superk<kt> superk("CATACAGAGACAGCAGCAGAGCA", 20);
  cr_assert(superk.str_value() == "CATACAGAGACAGCAGCAGAGCA");
  uchar buf[8] = {0x11, 0x73, 0x40, 0x34, 0xCD, 0x00, 0xC0, 0x0C};
  superk.set_superk(buf, 30, 20, true);
  cr_assert(superk.str_value() == "GAGCAGCACAAACGAGACACAAAAAAAGAG");
}

Test(superk, superk_get_kmer)
{
  Superk<kt> superk("GAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  cr_assert(superk.get_kmer(0, false).value() == 0xCD34407311);
  cr_assert(superk.get_kmer(1, false).value() == 0x34D101CC44);
}

Test(superk, superk_nb_kmers)
{
  Superk<kt> superk("GAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  cr_assert(superk.nb_kmers() == 11);
}

Test(minimizer, minimizer_build_from_kmer)
{
  Kmer<kt> kmer("GAGCAGCACAAACGAGACAC", true);
  Minimizer<kt> minim(&kmer, 10, false);
  cr_assert(minim.value() == 0x1CC4);
  cr_assert(minim.str_value() == "AAACGAGACA"); // not valid minimizer
  Minimizer<kt> minim2(&kmer, 10, true); // check validity
  cr_assert(minim2.value() == 0x7311);
  cr_assert(minim2.str_value() == "AACGAGACAC");
}

Test(minimizer, minimizer_build_from_superk)
{
  Superk<kt> superk("GAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  Minimizer<kt> minim(&superk, 10, false);
  cr_assert(minim.value() == 0x1CC4);
  cr_assert(minim.str_value() == "AAACGAGACA");
}

Test(minimizer, minimizer_set_from_kmer)
{
  Kmer<kt> kmer("GAGCAGCACAAACGAGACAC", true);
  Minimizer<kt> minim(10);
  minim.set_kmer(&kmer, 10, true);
  cr_assert(minim.value() == 0x7311);
}

Test(minimizer, minimizer_set_from_superk)
{
  Superk<kt> superk("GAGCAGCACAAACGAGACACAAAAAAAGAG", 20);
  Minimizer<kt> minim(10);
  minim.set_superk(&superk, 10, true);
  cr_assert(minim.value() == 0x7311);
}

Test(minimizer, minimizer_operator)
{
  Kmer<kt> kmer("GAGCAGCACAAACGAGACAC", true);
  Minimizer<kt> minim(&kmer, 10, false);

  cr_assert(minim < 0xFFFF);
  cr_assert(minim > 0xFF);
  cr_assert(minim != 0xFF);
  cr_assert(minim == 0x1CC4);

  cr_assert(minim < "AGCGTACGAA");
  cr_assert(minim > "AAAAACGTAA");
  cr_assert(minim != "ACGTACGGAA");
  cr_assert(minim == "AAACGAGACA");
}

Test(minimizer, minimizer_default)
{
  Kmer<kt> kmer("AAAAAAAAAAAAAAAAAAAAAAA", true);
  Minimizer<kt> minim(&kmer, 10, true);
  cr_assert(minim.value() == DEFAULT_MINIMIZER);

  minim.set_default(0x1234);
  cr_assert(minim.value() == 0x1234);

  minim.set_default("ACGTACGTAA");
  cr_assert(minim.value() == 0x01E1E0);

  minim.set_default();
  cr_assert(minim.value() == DEFAULT_MINIMIZER);
}