#include <gtest/gtest.h>
#include <kmtricks/io/fof.hpp>

using namespace km;

TEST(fof, fof)
{
  Fof fof("./data/fof.txt");

  auto it = fof.begin();
  EXPECT_EQ(std::get<0>(*it), "D1");
  EXPECT_EQ(std::get<2>(*it), 20);
  EXPECT_EQ(std::get<1>(*it)[0], "/path/to/D1.fasta");
  EXPECT_EQ(std::get<1>(*it)[1], "/path/to/D1.fasta");
  it++;
  EXPECT_EQ(std::get<0>(*it), "D2");
  EXPECT_EQ(std::get<2>(*it), 20);
  EXPECT_EQ(std::get<1>(*it)[0], "/path/to/D2.fasta");
  EXPECT_EQ(std::get<1>(*it)[1], "/path/to/D2.fasta");
}