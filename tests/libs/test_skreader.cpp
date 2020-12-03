#include <gtest/gtest.h>
#define _KM_LIB_INCLUDE_
#include <kmtricks/skreader.hpp>
#include <kmtricks/sequences.hpp>

using namespace km;
typedef uint64_t kt;

TEST(skstorage, skstorage_build)
{
  string path = "./skreader_data/sk_part";
  string prefix = "superKparts.";
  SuperkStorage s(path, prefix, 4);

  EXPECT_EQ(s.nb_files(), 4);

  for (int i=0; i<s.nb_files(); i++)
  {
    EXPECT_TRUE(s._parts[i].is_open());
  }
}

TEST(skstorage, skstorage_read_block)
{
  string path = "./skreader_data/sk_part";
  string prefix = "superKparts.";
  SuperkStorage s(path, prefix, 4);

  uchar *buffer = 0;
  uint  size = 0;
  uint  nb_bytes = 0;

  int r = s.read_block(&buffer, &size, &nb_bytes, 0);
  EXPECT_EQ(r, 72);
  EXPECT_EQ(*buffer, 4);

  Superk<kt> superk(++buffer, 23, 20, true);
  EXPECT_EQ(superk.str_value(), "CATACAGAGACAGCAGCAGAGCA");
}

TEST(skreader, skreader_build)
{
  string path = "./skreader_data/sk_part";
  string prefix = "superKparts.";
  SuperkStorage s(path, prefix, 4);
  size_t kmersize = 20;
  SuperkReader<kt> reader(&s, kmersize);
  Superk<kt> superk(kmersize);
}

TEST(skreader, skreader_next)
{
  string path = "./skreader_data/sk_part";
  string prefix = "superKparts.";
  SuperkStorage s(path, prefix, 4);
  size_t kmersize = 20;
  SuperkReader<kt> reader(&s, kmersize);
  Superk<kt> superk(kmersize);
  reader.next_superk(0, &superk);
  EXPECT_EQ(superk.str_value(), "CATACAGAGACAGCAGCAGAGCA");
  reader.next_superk(0, &superk);
  EXPECT_EQ(superk.str_value(), "GAGCAGCACAAACGAGACACAAAAAAAGAG");
}