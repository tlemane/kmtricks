#include <gtest/gtest.h>
#include <kmtricks/merge.hpp>


TEST(merge, hash_merge)
{
  std::vector<std::vector<std::string>> paths;
  for (size_t i=0; i<4; i++)
  {
    std::vector<std::string> p = {
      "./data/partitions/hashes/partition_" + std::to_string(i) + "/D1.hash",
      "./data/partitions/hashes/partition_" + std::to_string(i) + "/D2.hash",
    };
    paths.push_back(p);
  }
  std::vector<uint32_t> a {1, 1};
  {
    int count = 0;
    km::HashMerger<255, 32768, km::HashReader<255>> m(paths[0], a, 1, 1);
    while (m.next()) { count++; }
    EXPECT_EQ(count, 57);
  }
  {
    int count = 0;
    km::HashMerger<255, 32768, km::HashReader<255>> m(paths[1], a, 1, 1);
    while (m.next()) { count++; }
    EXPECT_EQ(count, 67);
  }
  {
    int count = 0;
    km::HashMerger<255, 32768, km::HashReader<255>> m(paths[2], a, 1, 1);
    while (m.next()) { count++; }
    EXPECT_EQ(count, 70);
  }
  {
    int count = 0;
    km::HashMerger<255, 32768, km::HashReader<255>> m(paths[3], a, 1, 1);
    while (m.next()) { count++; }
    EXPECT_EQ(count, 82);
  }
}

TEST(merge, kmer_merge)
{
  std::vector<std::vector<std::string>> paths;
  for (size_t i=0; i<4; i++)
  {
    std::vector<std::string> p = {
      "./data/partitions/kmers/partition_" + std::to_string(i) + "/D1.kmer",
      "./data/partitions/kmers/partition_" + std::to_string(i) + "/D2.kmer",
    };
    paths.push_back(p);
  }
  std::vector<uint32_t> a {1, 1};
  {
    int count = 0;
    km::KmerMerger<32, std::numeric_limits<uint32_t>::max()> m(paths[0], a, 31, 1, 1);
    while (m.next()) { count++; }
    EXPECT_EQ(count, 57);
  }
  {
    int count = 0;
    km::KmerMerger<32, std::numeric_limits<uint32_t>::max()> m(paths[1], a, 31, 1, 1);
    while (m.next()) { count++; }
    EXPECT_EQ(count, 67);
  }
  {
    int count = 0;
    km::KmerMerger<32, std::numeric_limits<uint32_t>::max()> m(paths[2], a, 31, 1, 1);
    while (m.next()) { count++; }
    EXPECT_EQ(count, 70);
  }
  {
    int count = 0;
    km::KmerMerger<32, std::numeric_limits<uint32_t>::max()> m(paths[3], a, 31, 1, 1);
    while (m.next()) { count++; }
    EXPECT_EQ(count, 82);
  }
}