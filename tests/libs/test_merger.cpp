#include <gtest/gtest.h>
#define _KM_LIB_INCLUDE_
#include <kmtricks/merger.hpp>

using namespace km;
typedef uint64_t kt;
typedef uint8_t  kc;

TEST(merger, merger_build)
{
  string path = "./merge_data/fof.txt";
  Merger<kt, kc> m(path, 1, 1, 0, true);
  EXPECT_EQ(m.nb_files, 2);
}

TEST(merger, merger_next)
{
  string path = "./merge_data/fof.txt";
  Merger<kt, kc> m(path, 1, 2, 0, true);

  kt values[5] = {0, 134, 234, 300, 302};
  kc count1[5] = {1, 31, 1, 8, 12};
  kc count2[5] = {3, 6, 100, 2, 1};

  int i = 0;
  while(!m.end)
  {
    m.next();
    if (m.keep)
    {
      EXPECT_EQ(m.m_khash, values[i]);
      EXPECT_EQ(m.counts[0], count1[i]);
      EXPECT_EQ(m.counts[1], count2[i]);
      i++;
    }
  }
}

TEST(merger, merger_get_kmer)
{
  string path = "./merge_data/fof.txt";
  Merger<kt, kc> m(path, 1, 2, 0, true);

  kt values[5] = {0, 134, 234, 300, 302};

  int i = 0;
  while(!m.end)
  {
    m.next();
    if (m.keep)
    {
      Kmer<kt> k(values[i], 32, false);
      Kmer<kt> k2 = m.get_kmer(32);
      EXPECT_EQ(k, k2);
      i++;
    }
  }
}