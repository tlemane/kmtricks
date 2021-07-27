#include <gtest/gtest.h>
#define KMER_LIST 32, 64, 96
#define KMER_N 3

#include <kmtricks/loop_executor.hpp>

template<size_t M>
struct TestFunctor
{
  void operator()(int i, size_t& value)
  {
    value = M;
  }
};

TEST(kmer_selector, kmer_size_selector)
{
  size_t value = 0;
  km::const_loop_executor<0, KMER_N>::exec<TestFunctor>(30, 1, value);
  EXPECT_EQ(value, 32);
  km::const_loop_executor<0, KMER_N>::exec<TestFunctor>(60, 42, value);
  EXPECT_EQ(value, 64);
  km::const_loop_executor<0, KMER_N>::exec<TestFunctor>(90, 42, value);
  EXPECT_EQ(value, 96);
}