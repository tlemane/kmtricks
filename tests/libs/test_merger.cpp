#include <gtest/gtest.h>
#include <string>
#include <fstream>
#include <vector>
#include <cstdio>
#define _KM_LIB_INCLUDE_
#include <kmtricks/merger.hpp>
#include <kmtricks/io.hpp>

using namespace km;
using namespace std;
typedef uint64_t kt;
typedef uint8_t  kc;


class MTest : public ::testing::Test
{
protected:
  MTest ()
  {
    Kmer<kt> k(false);
    pathInt = {"./merge_data/toy1.txt", "./merge_data/toy2.txt"};
    pathKmer = {"./merge_data/toy1_kmer.txt", "./merge_data/toy2_kmer.txt"};
    pathBin = {"./merge_data/toy1_kmer.txt.bin", "./merge_data/toy2_kmer.txt.bin"};

    for (size_t i=0; i<pathInt.size(); i++)
    {
      ifstream f1(pathInt[i], ios::in);
      ofstream f1_txt(pathKmer[i], ios::out);
      KmerFile<OUT, kt, kc> f1bin(
        pathBin[i], 0, 0, 10, 0, 0);
      string line; 
      while (getline(f1, line))
      {
        auto v = split(line, ' ');
        uint64_t ku = stoi(v[0]);
        uint8_t cu = stoi(v[1]);
        k.set_kmer(ku, 10);
        f1_txt << k.str_value() << " " << to_string(cu) << "\n";
        f1bin.write(k, cu);
      }
    }
  }
  ~MTest ()
  {
  }

  vector<string> pathInt;
  vector<string> pathKmer;
  vector<string> pathBin;
};

TEST_F(MTest, merger_build)
{
  string path = "./merge_data/fof_kmer.txt";
  Merger<kt, kc, KmerFile<IN, kt, kc>> m(path, 1, 1, 0, true);
  EXPECT_EQ(m.nb_files, 2);
}

TEST_F(MTest, merger_next)
{
  string path = "./merge_data/fof_kmer.txt";
  Merger<kt, kc, KmerFile<IN, kt, kc>> m(path, 1, 2, 0, true);

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

TEST_F(MTest, merger_get_kmer)
{
  string path = "./merge_data/fof_kmer.txt";
  Merger<kt, kc, KmerFile<IN, kt, kc>> m(path, 1, 2, 0, true);

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