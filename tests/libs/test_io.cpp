#include <random>
#include <vector>
#include <unordered_map>
#include <cstdio>
#include <gtest/gtest.h>
#define _KM_LIB_INCLUDE_
#include <kmtricks/io.hpp>

using ktype = uint64_t;
using cntype = uint16_t;

using namespace km;
using namespace std;
auto file_id        = 1;
auto partition_id   = 1;
auto partition_size = 1024;
auto cols           = 120;
auto kmer_size      = 10;
auto c0             = 0;
auto c1             = 1;
auto N              = 1024;
auto C              = cols;

vector<cntype> get_rand(int n, uniform_int_distribution<cntype> distC, mt19937 mt)
{
  vector<cntype> v(n);
  for (int i=0; i<n; i++)
    v[i] = distC(mt);
  return v;
}
class IoTest : public ::testing::Test
{
protected:
  IoTest () {
    mt19937 mt(rd());
    for (int i=0; i<N; i++)
    {
      kmap.insert(
        {Kmer<ktype>(distK(mt), kmer_size, false).value(), get_rand(C, distC, mt)}
      );
      vector<char> bits(NBYTE(C));
      for (auto& byte : bits) byte = dist8(mt);
      kbit.insert(
        {Kmer<ktype>(distK(mt), kmer_size, false).value(), bits}
      );
      bfs.push_back(bits);
    }
    for (int i=0; i<1024; i++)
      for (int j=0; j<120; j++ )
        mat.set_bit(i, 0, 1);
    mat.set_bit(1, 0, 0);
  }

  void SetUp() override {
    
  }

  void TearDown() override {
 
 }

  ~IoTest()
  {
    std::remove("KmerFile");
    std::remove("CountMatrixASCII");
    std::remove("CountMatrixBIN");
    std::remove("PAMatrix");
    std::remove("BitMatrixFileBF");
    std::remove("MatrixFileBIT");
    std::remove("BitVectorFile");
    std::remove("HistFile");
  }

  unordered_map<ktype, vector<cntype>> kmap;
  unordered_map<ktype, vector<char>> kbit;
  vector<vector<char>> bfs;
  random_device rd;
  uniform_int_distribution<ktype> distK;
  uniform_int_distribution<cntype> distC;
  uniform_int_distribution<char> dist8;
  BitMatrix mat {BitMatrix(1024, 15, true)};
};


TEST(io, KmerFileFile)
{
  ktype kmer = 1024;
  cntype count = 512;
  {
    KmerFile<OUT, ktype, cntype> f("KmerFile", file_id, partition_id, kmer_size, 0, c0);
    f.write(kmer, count);
  }

  {
    std::ifstream fr("KmerFile", std::ios::binary | std::ios::in);
    ktype kmer_read;
    cntype count_read;
    fr.seekg(sizeof(kheader_t));
    fr.read(reinterpret_cast<char*>(&kmer_read), sizeof(ktype));
    fr.read(reinterpret_cast<char*>(&count_read), sizeof(cntype));
    EXPECT_EQ(kmer_read, kmer);
    EXPECT_EQ(count_read, count);
  }
  {
    KmerFile<IN, ktype, cntype> fr("KmerFile");
    ktype kmer_read;
    cntype count_read;
    fr.read(kmer_read, count_read);
    EXPECT_EQ(kmer_read, kmer);
    EXPECT_EQ(count_read, count);
  }
}

TEST_F(IoTest, CountMatrixFile)
{
  {
    CountMatrixFile<OUT, ktype, cntype, matrix_format_t::ASCII> asciiM(
      "CountMatrixASCII", partition_id, cols, kmer_size, 0, c0
    );
    CountMatrixFile<OUT, ktype, cntype, matrix_format_t::BIN> binM(
      "CountMatrixBIN", partition_id, cols, kmer_size, 0, c0
    );
    Kmer<ktype> k(false);
    for (auto& kc : kmap)
    {
      k.set_kmer(kc.first, kmer_size);
      asciiM.write(k, kc.second);
      binM.write(k, kc.second);
    }
  }
  {
    CountMatrixFile<IN, ktype, cntype, matrix_format_t::ASCII> asciiM(
      "CountMatrixASCII"
    );
    CountMatrixFile<IN, ktype, cntype, matrix_format_t::BIN> binM(
      "CountMatrixBIN"
    );
    Kmer<ktype> k(false);
    vector<cntype> counts(cols);
    while(asciiM.read(k, counts))
    {
      EXPECT_TRUE(kmap.count(k.value()));
      for (size_t i=0; i<counts.size(); i++)
        EXPECT_EQ(counts[i], kmap.at(k.value())[i]);
    }
    while(binM.read(k, counts))
    {
      EXPECT_TRUE(kmap.count(k.value()));
      for (size_t i=0; i<counts.size(); i++)
        EXPECT_EQ(counts[i], kmap.at(k.value())[i]);
    }
  }
}

TEST_F(IoTest, PAMatrixFile)
{
  {
    PAMatrixFile<OUT, ktype> pam(
      "PAMatrix", partition_id, cols, kmer_size, 0, c0
    );
    Kmer<ktype> k(false);
    for (auto& kc : kbit)
    {
      k.set_kmer(kc.first, kmer_size);
      pam.write(k, kc.second);
    }
  }
  {
    PAMatrixFile<IN, ktype> pam(
      "PAMatrix"
    );
    Kmer<ktype> k(false);
    vector<char> counts(pam.infos()->size_in_bytes);
    while(pam.read(k, counts))
    {
      EXPECT_TRUE(kbit.count(k.value()));
      for (size_t i=0; i<counts.size(); i++)
        EXPECT_EQ(counts[i], kbit.at(k.value())[i]);
    }
  }
}

TEST_F(IoTest, matrix_file_bf)
{
  {
    BitMatrixFile<ostream, matrix_format_t::BF> matrix(
      "BitMatrixFileBF", partition_id, 1024, 120, c0);
    for (auto& v: bfs)
    {
      EXPECT_FALSE(matrix.is_consistent());
      EXPECT_TRUE(matrix.write(v));
    }
    EXPECT_TRUE(matrix.is_consistent());
    EXPECT_FALSE(matrix.write(bfs[0])); // add warning for this
  }
  {
    BitMatrixFile<istream, matrix_format_t::BF> matrix(
      "BitMatrixFileBF");
    vector<char> bits(matrix.infos()->nb_cols/8);
    int i = 0;
    while (matrix.read(bits))
    {
      auto& v = bfs[i];
      EXPECT_TRUE(equal(v.begin(), v.end(), bits.begin()));
      i++;
    }
  }
}

TEST_F(IoTest, BitMatrixFile)
{
  {
    BitMatrix *trp = mat.transpose();
    auto i = trp->get_nb_lines();
    auto j = trp->get_nb_cols();
    BitMatrixFile<ostream, matrix_format_t::BIT> MF(
      "MatrixFileBIT", partition_id, i, j, c0);
    MF.dump(*trp);
    delete trp;
  }
  {
    BitMatrixFile<istream, matrix_format_t::BIT> MF(
      "MatrixFileBIT"
    );
    auto i = MF.infos()->nb_rows_use;
    auto j = MF.infos()->nb_cols_use/8;
    BitMatrix mat_read(i, j, true);
    MF.load(mat_read);
    BitMatrix *retrp = mat_read.transpose();
    EXPECT_TRUE(equal(mat.matrix, mat.matrix+mat.get_size_in_byte(), retrp->matrix));
    delete retrp;
  }
}

TEST_F(IoTest, BitVectorFile)
{
  {
    BitVectorFile<ostream> bvfile(
      "BitVectorFile", file_id, partition_id, cols, c0
    );
    bvfile.write(bfs[0]);
  }
  {
    BitVectorFile<istream> bvfile(
      "BitVectorFile"
    );
    vector<char> v = std::move(bvfile.read());
    EXPECT_TRUE(equal(bfs[0].begin(), bfs[0].end(), v.begin()));

    pair<uint64_t, uint64_t> w = bvfile.get_window();
    EXPECT_EQ(w.first, 120);
    EXPECT_EQ(w.second, 239);
  }
}

TEST_F(IoTest, HistFile)
{
  vector<uint64_t> v {1, 1, 3, 9, 1, 2, 2, 2, 9, 5};
  vector<uint64_t> r {3, 3, 1, 0, 1, 0, 0, 0, 2, 0};
  vector<uint64_t> rn {3, 6, 3, 0, 5, 0, 0, 0, 18, 0};
  {
    KHist hist(0, 20, 1, 10);

    for (auto& c : v)
    {
      hist.inc(c);
    }
    HistFile<OUT> hf(hist, "HistFile");
  }
  {
    HistFile<IN> hf("HistFile");
    KHist hist = hf.read();
    hist.print_histu();
    EXPECT_EQ(hist.lower, 1);
    EXPECT_EQ(hist.upper, 10);
    EXPECT_EQ(hist.oob_un, 0);
    EXPECT_EQ(hist.oob_ln, 0);
    EXPECT_EQ(hist.oob_lu, 0);
    EXPECT_EQ(hist.oob_uu, 0);
    for (int i=0; i<v.size(); i++)
    {
      EXPECT_EQ(r[i], hist.hist_u[i]);
      EXPECT_EQ(rn[i], hist.hist_n[i]);
    }
  }
}