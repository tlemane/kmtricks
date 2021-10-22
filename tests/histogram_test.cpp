#include <gtest/gtest.h>
#include <kmtricks/histogram.hpp>
#include <kmtricks/io/hist_file.hpp>

using namespace km;

TEST(histogram, histogram)
{
  std::vector<uint64_t> v {1, 1, 3, 9, 1, 2, 2, 2, 9, 5};
  std::vector<uint64_t> r {3, 3, 1, 0, 1, 0, 0, 0, 2, 0};
  std::vector<uint64_t> rn {3, 6, 3, 0, 5, 0, 0, 0, 18, 0};
  {
    KHist hist(0, 20, 1, 10);

    for (auto& c : v)
    {
      hist.inc(c);
    }
    HistWriter<8192> hw("./tests_tmp/h.hist", hist, false);
  }
  {
    HistReader<8192> hr("./tests_tmp/h.hist");
    hist_t hist = hr.get();
    EXPECT_EQ(hist->lower(), 1);
    EXPECT_EQ(hist->upper(), 10);
    EXPECT_EQ(hist->oob_upper_total(), 0);
    EXPECT_EQ(hist->oob_lower_total(), 0);
    EXPECT_EQ(hist->oob_lower_unique(), 0);
    EXPECT_EQ(hist->oob_upper_unique(), 0);
    for (int i=0; i<v.size(); i++)
    {
      EXPECT_EQ(r[i], hist->get_vec(KHistType::UNIQUE)[i]);
      EXPECT_EQ(rn[i], hist->get_vec(KHistType::TOTAL)[i]);
    }
  }
}

TEST(histogram, clones)
{
  std::vector<uint64_t> v {1, 1, 3, 9, 1};
  std::vector<uint64_t> v2 {2, 2, 2, 9, 5};

  std::vector<uint64_t> r {3, 3, 1, 0, 1, 0, 0, 0, 2, 0};
  std::vector<uint64_t> rn {3, 6, 3, 0, 5, 0, 0, 0, 18, 0};

  hist_t hist = std::make_shared<KHist>(0, 20, 1, 10);

  hist_t c1 = hist->clone();
  hist_t c2 = hist->clone();

  for (auto& c : v)
    c1->inc(c);
  for (auto& c : v2)
    c2->inc(c);

  hist->merge_clones();

  EXPECT_EQ(hist->lower(), 1);
  EXPECT_EQ(hist->upper(), 10);
  EXPECT_EQ(hist->oob_upper_total(), 0);
  EXPECT_EQ(hist->oob_lower_total(), 0);
  EXPECT_EQ(hist->oob_lower_unique(), 0);
  EXPECT_EQ(hist->oob_upper_unique(), 0);
  for (int i=0; i<v.size(); i++)
  {
    EXPECT_EQ(r[i], hist->get_vec(KHistType::UNIQUE)[i]);
    EXPECT_EQ(rn[i], hist->get_vec(KHistType::TOTAL)[i]);
  }
}
