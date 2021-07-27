#include <gtest/gtest.h>
#include <kmtricks/io/pa_matrix_file.hpp>
#include <kmtricks/utils.hpp>

using namespace km;

TEST(matrix_file, PAMatrixWriter)
{
  PAMatrixWriter pw("tests_tmp/p1.matrix", 21, 20, 1, 2, false);
  EXPECT_EQ(pw.infos().kmer_size, 21);
  EXPECT_EQ(pw.infos().kmer_slots, (21+31)/32);
  EXPECT_EQ(pw.infos().bits, 20);
  EXPECT_EQ(pw.infos().bytes, NBYTES(20));
  EXPECT_EQ(pw.infos().id, 1);
  EXPECT_EQ(pw.infos().partition, 2);
  EXPECT_EQ(pw.infos().compressed, false);
}

TEST(matrix_file, PAMatrixReader)
{
  PAMatrixReader pr("tests_tmp/p1.matrix");
  EXPECT_EQ(pr.infos().kmer_size, 21);
  EXPECT_EQ(pr.infos().kmer_slots, (21+31)/32);
  EXPECT_EQ(pr.infos().bits, 20);
  EXPECT_EQ(pr.infos().bytes, NBYTES(20));
  EXPECT_EQ(pr.infos().id, 1);
  EXPECT_EQ(pr.infos().partition, 2);
  EXPECT_EQ(pr.infos().compressed, false);
}

TEST(matrix_file, PAMatrixWriteRead)
{
  std::vector<std::string> str_kmers(10000);
  std::vector<std::vector<uint8_t>> counts(10000, std::vector<uint8_t>(NBYTES(20)));
  {
    PAMatrixWriter pw("tests_tmp/p2.matrix", 21, 20, 1, 2, false);
    PAMatrixWriter pw2("tests_tmp/p2.matrix.lz4", 21, 20, 1, 2, true);
    for (size_t i=0; i<str_kmers.size(); i++)
    {
      str_kmers[i] = random_dna_seq(21);
      counts[i] = random_count_vector<uint8_t>(NBYTES(20));
      Kmer<32> kmer(str_kmers[i]);
      pw.write<32>(kmer, counts[i]);
      pw2.write<32>(kmer, counts[i]);
    }
  }
  {
    PAMatrixReader pw("tests_tmp/p2.matrix");
    PAMatrixReader pw2("tests_tmp/p2.matrix.lz4");
    Kmer<32> kmer; kmer.set_k(pw.infos().kmer_size);
    std::vector<uint8_t> c(pw.infos().bytes);

    for (size_t i=0; i<str_kmers.size(); i++)
    {
      pw.read<32>(kmer, c);
      EXPECT_EQ(kmer.to_string(), str_kmers[i]);
      EXPECT_TRUE(std::equal(c.begin(), c.end(), counts[i].begin()));
      pw2.read<32>(kmer, c);
      EXPECT_EQ(kmer.to_string(), str_kmers[i]);
      EXPECT_TRUE(std::equal(c.begin(), c.end(), counts[i].begin()));
    }
  }
  {
    std::ofstream out("tests_tmp/p2.matrix.csv");
    PAMatrixReader("tests_tmp/p2.matrix").write_as_text<32>(out);
  }
}