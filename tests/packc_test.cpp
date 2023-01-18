#include <gtest/gtest.h>

#include <kmtricks/packc.hpp>

TEST(packc, byte_count_pack)
{
  EXPECT_EQ(km::byte_count_pack(10, 5), 7);
  EXPECT_EQ(km::byte_count_pack(10, 10), 13);
  EXPECT_EQ(km::byte_count_pack(10, 12), 15);
  EXPECT_EQ(km::byte_count_pack(0, 0), 0);
  EXPECT_EQ(km::byte_count_pack(0, 1), 0);
  EXPECT_EQ(km::byte_count_pack(100, 0), 0);
  EXPECT_EQ(km::byte_count_pack(10, 2), 3);
  EXPECT_EQ(km::byte_count_pack(8, 2), 2);
  EXPECT_EQ(km::byte_count_pack(2, 8), 2);
  EXPECT_EQ(km::byte_count_pack(3, 8), 3);
}

TEST(packc, to_n_b)
{
  // base case
  EXPECT_EQ(km::to_n_b(0, 100), 0);
  EXPECT_EQ(km::to_n_b(1, 100), 1);
  EXPECT_EQ(km::to_n_b(2, 100), 2);
  EXPECT_EQ(km::to_n_b(3, 100), 2);
  EXPECT_EQ(km::to_n_b(4, 100), 3);
  EXPECT_EQ(km::to_n_b(5, 100), 3);
  EXPECT_EQ(km::to_n_b(6, 100), 3);
  EXPECT_EQ(km::to_n_b(7, 100), 3);
  EXPECT_EQ(km::to_n_b(8, 100), 4);
  EXPECT_EQ(km::to_n_b(9, 100), 4);
  EXPECT_EQ(km::to_n_b(16, 100), 5);
  EXPECT_EQ(km::to_n_b(32, 100), 6);
  EXPECT_EQ(km::to_n_b(32, 2), 3);  // caped
  EXPECT_EQ(km::to_n_b(32767, 100), 15);
  EXPECT_EQ(km::to_n_b(32768, 100), 16);
  EXPECT_EQ(km::to_n_b(32769, 100), 16);
  EXPECT_EQ(km::to_n_b(32767, 3), 7);  // caped
  EXPECT_EQ(km::to_n_b(32768, 3), 7);  // caped
  EXPECT_EQ(km::to_n_b(32769, 3), 7);  // caped
}