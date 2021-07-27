#include <gtest/gtest.h>
#include <kmtricks/io/hash_file.hpp>
#include <kmtricks/utils.hpp>

using namespace km;

TEST(hash_file, HashWriter)
{
  HashWriter<255> kw("tests_tmp/h1.hash", 1, 1, 2, false);
  EXPECT_EQ(kw.infos().count_slots, 1);
  EXPECT_EQ(kw.infos().id, 1);
  EXPECT_EQ(kw.infos().partition, 2);
  EXPECT_EQ(kw.infos().compressed, false);
}

TEST(hash_file, HashReader)
{
  HashReader<255> kr("tests_tmp/h1.hash");
  EXPECT_EQ(kr.infos().count_slots, 1);
  EXPECT_EQ(kr.infos().id, 1);
  EXPECT_EQ(kr.infos().partition, 2);
  EXPECT_EQ(kr.infos().compressed, false);
}

TEST(hash_file, HashWriteRead)
{
  {
    HashWriter<255> kw("tests_tmp/h2.hash", 1, 1, 2, false);
    HashWriter<255> kw2("tests_tmp/h2.hash.lz4", 1, 1, 2, true);
    for (uint64_t i=0; i<10000; i++)
    {
      kw.write(i, 42);
      kw2.write(i, 42);
    }
  }
  {
    HashReader<255> kr("tests_tmp/h2.hash");
    HashReader<255> kr2("tests_tmp/h2.hash.lz4");
    uint64_t hash;
    uint8_t c = 0;
    for (uint64_t i=0; i<10000; i++)
    {
      kr.read(hash, c);
      EXPECT_EQ(hash, i);
      EXPECT_EQ(c, 42);
      kr2.read(hash, c);
      EXPECT_EQ(hash, i);
      EXPECT_EQ(c, 42);
    }
  }
}
