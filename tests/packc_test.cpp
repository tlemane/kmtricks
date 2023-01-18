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
  EXPECT_EQ(km::to_n_b(0, 8), 0);
  EXPECT_EQ(km::to_n_b(1, 8), 1);
  EXPECT_EQ(km::to_n_b(2, 8), 2);
  EXPECT_EQ(km::to_n_b(3, 8), 2);
  EXPECT_EQ(km::to_n_b(4, 8), 3);
  EXPECT_EQ(km::to_n_b(5, 8), 3);
  EXPECT_EQ(km::to_n_b(6, 8), 3);
  EXPECT_EQ(km::to_n_b(7, 8), 3);
  EXPECT_EQ(km::to_n_b(8, 8), 4);
  EXPECT_EQ(km::to_n_b(9, 8), 4);
  EXPECT_EQ(km::to_n_b(16, 8), 5);
  EXPECT_EQ(km::to_n_b(32, 8), 6);
  EXPECT_EQ(km::to_n_b(32, 2), 3);  // caped
  EXPECT_EQ(km::to_n_b(32767, 8), 15);
  EXPECT_EQ(km::to_n_b(32768, 8), 16);
  EXPECT_EQ(km::to_n_b(32769, 8), 16);
  EXPECT_EQ(km::to_n_b(32767, 3), 7);  // caped
  EXPECT_EQ(km::to_n_b(32768, 3), 7);  // caped
  EXPECT_EQ(km::to_n_b(32769, 3), 7);  // caped
}