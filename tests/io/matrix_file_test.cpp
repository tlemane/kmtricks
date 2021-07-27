#include <gtest/gtest.h>
#include <kmtricks/io/matrix_file.hpp>
#include <kmtricks/utils.hpp>

using namespace km;

TEST(matrix_file, MatrixWriter)
{
  MatrixWriter mw("tests_tmp/m1.matrix", 21, 1, 10, 1, 2, false);
  EXPECT_EQ(mw.infos().kmer_size, 21);
  EXPECT_EQ(mw.infos().kmer_slots, (21+31)/32);
  EXPECT_EQ(mw.infos().count_slots, 1);
  EXPECT_EQ(mw.infos().nb_counts, 10);
  EXPECT_EQ(mw.infos().id, 1);
  EXPECT_EQ(mw.infos().partition, 2);
  EXPECT_EQ(mw.infos().compressed, false);
}

TEST(matrix_file, MatrixReader)
{
  MatrixReader mr("tests_tmp/m1.matrix");
  EXPECT_EQ(mr.infos().kmer_size, 21);
  EXPECT_EQ(mr.infos().kmer_slots, (21+31)/32);
  EXPECT_EQ(mr.infos().count_slots, 1);
  EXPECT_EQ(mr.infos().nb_counts, 10);
  EXPECT_EQ(mr.infos().id, 1);
  EXPECT_EQ(mr.infos().partition, 2);
  EXPECT_EQ(mr.infos().compressed, false);
}

TEST(matrix_file, MatrixWriteRead)
{
  std::vector<std::string> str_kmers(10000);
  std::vector<std::vector<uint8_t>> counts(10000, std::vector<uint8_t>(50));
  {
    MatrixWriter mw("tests_tmp/m2.matrix", 21, 1, 50, 1, 2, false);
    MatrixWriter mw2("tests_tmp/m2.matrix.lz4", 21, 1, 50, 1, 2, true);
    for (size_t i=0; i<str_kmers.size(); i++)
    {
      str_kmers[i] = random_dna_seq(21);
      counts[i] = random_count_vector<uint8_t>(50);
      Kmer<32> kmer(str_kmers[i]);
      mw.write<32, 255>(kmer, counts[i]);
      mw2.write<32, 255>(kmer, counts[i]);
    }
  }
  {
    MatrixReader rw("tests_tmp/m2.matrix");
    MatrixReader rw2("tests_tmp/m2.matrix.lz4");
    Kmer<32> kmer; kmer.set_k(rw.infos().kmer_size);
    std::vector<uint8_t> c(rw.infos().nb_counts);

    for (size_t i=0; i<str_kmers.size(); i++)
    {
      rw.read<32, 255>(kmer, c);
      EXPECT_EQ(kmer.to_string(), str_kmers[i]);
      EXPECT_TRUE(std::equal(c.begin(), c.end(), counts[i].begin()));
      rw2.read<32, 255>(kmer, c);
      EXPECT_EQ(kmer.to_string(), str_kmers[i]);
      EXPECT_TRUE(std::equal(c.begin(), c.end(), counts[i].begin()));
    }
  }
  {
    std::ofstream out("tests_tmp/m2.matrix.csv");
    MatrixReader("tests_tmp/m2.matrix").write_as_text<32, 255>(out);
  }
}

TEST(matrix_file, MatrixHashWriter)
{
  MatrixHashWriter mw("tests_tmp/m1.hash_matrix", 1, 10, 1, 2, false);
  EXPECT_EQ(mw.infos().count_slots, 1);
  EXPECT_EQ(mw.infos().nb_counts, 10);
  EXPECT_EQ(mw.infos().id, 1);
  EXPECT_EQ(mw.infos().partition, 2);
  EXPECT_EQ(mw.infos().compressed, false);
}

TEST(matrix_file, MatrixHashReader)
{
  MatrixHashReader mr("tests_tmp/m1.hash_matrix");
  EXPECT_EQ(mr.infos().count_slots, 1);
  EXPECT_EQ(mr.infos().nb_counts, 10);
  EXPECT_EQ(mr.infos().id, 1);
  EXPECT_EQ(mr.infos().partition, 2);
  EXPECT_EQ(mr.infos().compressed, false);
}

TEST(matrix_file, MatrixHashWriteRead)
{
  std::vector<std::vector<uint8_t>> counts(10000, std::vector<uint8_t>(50));
  {
    MatrixHashWriter mw("tests_tmp/m2.hash_matrix", 1, 50, 1, 2, false);
    MatrixHashWriter mw2("tests_tmp/m2.hash_matrix.lz4", 1, 50, 1, 2, true);
    for (uint64_t i=0; i<10000; i++)
    {
      mw.write<255>(i, counts[i]);
      mw2.write<255>(i, counts[i]);
    }
  }
  {
    MatrixHashReader rw("tests_tmp/m2.hash_matrix");
    MatrixHashReader rw2("tests_tmp/m2.hash_matrix.lz4");
    std::vector<uint8_t> c(rw.infos().nb_counts);
    uint64_t hash;
    for (uint64_t i=0; i<10000; i++)
    {
      rw.read<255>(hash, c);
      EXPECT_EQ(hash, i);
      EXPECT_TRUE(std::equal(c.begin(), c.end(), counts[i].begin()));
      rw2.read<255>(hash, c);
      EXPECT_EQ(hash, i);
      EXPECT_TRUE(std::equal(c.begin(), c.end(), counts[i].begin()));
    }
  }
}