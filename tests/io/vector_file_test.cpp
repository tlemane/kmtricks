#include <gtest/gtest.h>
#include <kmtricks/io/vector_file.hpp>
#include <kmtricks/io/vector_matrix_file.hpp>
#include <kmtricks/utils.hpp>

using namespace km;

TEST(vector_file, BitVecWriter)
{
  BitVectorWriter bvw("./tests_tmp/b1.vec", 10000, 0, 1, false);
  EXPECT_EQ(bvw.infos().bits, 10000);
  EXPECT_EQ(bvw.infos().partition, 1);
  EXPECT_EQ(bvw.infos().id, 0);
  EXPECT_EQ(bvw.infos().compressed, false);
}

TEST(vector_file, BitVecReader)
{
  BitVectorReader bvr("./tests_tmp/b1.vec");
  EXPECT_EQ(bvr.infos().bits, 10000);
  EXPECT_EQ(bvr.infos().partition, 1);
  EXPECT_EQ(bvr.infos().id, 0);
  EXPECT_EQ(bvr.infos().compressed, false);
}

TEST(vector_file, BitVecReadWrite)
{
  std::vector<uint8_t> bits(NBYTES(10000));
  std::fill(bits.begin(), bits.end(), 42);
  {
    BitVectorWriter bvw("./tests_tmp/b2.vec", 10000, 0, 1, false);
    BitVectorWriter bvw2("./tests_tmp/b2.vec.lz4", 10000, 0, 1, true);
    bvw.write(bits);
    bvw2.write(bits);
  }
  {
    std::vector<uint8_t> tmp(NBYTES(10000));
    BitVectorReader bvr("./tests_tmp/b2.vec");
    BitVectorReader bvr2("./tests_tmp/b2.vec.lz4");
    bvr.read(tmp);
    EXPECT_TRUE(std::equal(bits.begin(), bits.end(), tmp.begin()));
    bvr2.read(tmp);
    EXPECT_TRUE(std::equal(bits.begin(), bits.end(), tmp.begin()));
  }
}