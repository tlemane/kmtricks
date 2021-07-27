#include <gtest/gtest.h>
#include <kmtricks/gatb/count_processor.hpp>

typedef ::Kmer<32>::Type Type;

auto revc = [](char c){
  return km::NToB[c];
};

TEST(processor, hash_count_processor)
{
  {
    km::hw_t<255> hw = std::make_shared<km::HashWriter<255>>("./tests_tmp/h.hash", 1, 0, 0, true);

    km::HashCountProcessor<32, 255> p(20, 3, hw, nullptr);
    p.process(0, 42, 2);
    p.process(0, 84, 6);
  }
  {
    uint64_t hash = 0;
    uint8_t c = 0;
    km::HashReader<255> hr("./tests_tmp/h.hash");
    hr.read(hash, c);
    EXPECT_EQ(hash, 84);
    EXPECT_EQ(c, 6);
    bool n = hr.read(hash, c);
    EXPECT_FALSE(n);
  }
}

TEST(processor, hash_vec_processor)
{
  {
    km::bvw_t<8192> hw = std::make_shared<km::BitVectorWriter<8192>>("./tests_tmp/hvec.hash", 1000, 0, 0, true);

    km::HashVecProcessor<32> p(20, 3, hw, nullptr, 1000);
    p.process(0, 42, 2);
    p.process(0, 84, 6);
    p.finish();
  }
  {
    std::vector<uint8_t> bits(NBYTES(1000));
    km::BitVectorReader hr("./tests_tmp/hvec.hash");
    hr.read(bits);
    EXPECT_TRUE(BITCHECK(bits, 84));
    EXPECT_FALSE(BITCHECK(bits, 42));
  }
}

TEST(processor, kmer_count_processor)
{
  std::string k1 = km::random_dna_seq(20);
  std::string k2 = km::random_dna_seq(20);
  Type gk1 = Type::polynom(k1.c_str(), 20, revc);
  Type gk2 = Type::polynom(k2.c_str(), 20, revc);
  {
    km::kw_t<8192> kw =
      std::make_shared<km::KmerWriter<8192>>("./tests_tmp/k.kmer", 20, 1, 0, 0, true);

    km::KmerCountProcessor<32, 255> p(20, 3, kw, nullptr);

    p.process(0, gk1, 2);
    p.process(0, gk2, 6);
  }
  {
    km::Kmer<32> kmer; kmer.set_k(20);
    uint8_t c = 0;
    km::KmerReader kr("./tests_tmp/k.kmer");
    kr.read<32, 255>(kmer, c);
    EXPECT_EQ(kmer.to_string(), gk2.toString(20));
    EXPECT_EQ(c, 6);
    bool a = kr.read<32, 255>(kmer, c);
    EXPECT_FALSE(a);
  }
}
