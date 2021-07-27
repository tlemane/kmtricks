#include <gtest/gtest.h>
#include <time.h>
#include <kmtricks/bitmatrix.hpp>
#include <kmtricks/io/vector_matrix_file.hpp>

typedef uint64_t kt;
using namespace km;

TEST(BitMatrix, bitmatrix_build)
{
  uint8_t zero[32]; memset(zero, 0, 32);
  uint8_t one[32]; memset(one, 0xFF, 32);
  BitMatrix mat(16, 2, true, false); // 16 line in bits, 2 cols in bytes -> 16*16
  EXPECT_TRUE(!memcmp(mat.matrix, zero, 32));
  BitMatrix mat1(16, 2, true, true);
  EXPECT_TRUE(!memcmp(mat1.matrix, one, 32));
}

TEST(BitMatrix, bitmatrix_clear)
{
  uint8_t zero[32]; memset(zero, 0, 32);
  BitMatrix mat(16, 2, true, true);
  mat.clear();
  EXPECT_TRUE(!(memcmp(mat.matrix, zero, 32)));
}

TEST(BitMatrix, bitmatrix_set_bit)
{
  BitMatrix mat(16, 2, true); // 16 line in bits, 2 cols in bytes -> 16*16
  EXPECT_EQ(mat.get_bit(4, 6), 0);
  mat.set_bit(4, 6, true);
  EXPECT_EQ(mat.get_bit(4, 6), 1);
}

TEST(BitMatrix, bitmatrix_tog_bit)
{
  BitMatrix mat(16, 2, true); // 16 line in bits, 2 cols in bytes -> 16*16
  mat.set_bit(4, 6, true);
  EXPECT_EQ(mat.get_bit(4, 6), 1);
  mat.tog_bit(4, 6);
  EXPECT_EQ(mat.get_bit(4, 6), 0);
}

TEST(BitMatrix, bitmatrix_set_byte)
{
  BitMatrix mat(16, 2, true);
  mat.set_byte(1, 1, 0x80);
  EXPECT_TRUE(mat.get_bit(1,15));
}

TEST(BitMatrix, bitmatrix_tog_byte)
{
  BitMatrix mat(16, 2, true);
  mat.set_byte(1, 1, 0x80);
  mat.tog_byte(1, 1);
  EXPECT_EQ(mat.get_byte(1, 1), 0x7F);
  EXPECT_FALSE(mat.get_bit(1,15));
}

TEST(BitMatrix, bitmatrix_tranpose)
{
  srand(time(NULL));
  BitMatrix mat(16, 2, true);
  for (int n=0; n<20; n++)
  {
    int i = rand() % 16;
    int j = rand() % 16;
    mat.set_bit(i, j, true);
  }

  BitMatrix *trp = mat.transpose();
  BitMatrix *rev = trp->transpose();

  EXPECT_TRUE(!memcmp(mat.matrix, rev->matrix, 32));


  {
    VectorMatrixWriter vmw("./tests_tmp/1.bit_matrix", 0, 0, 0, 0, 1, false);
    vmw.dump(mat);
  }
  {
    BitMatrix matl(16, 2, true);
    VectorMatrixReader vmw("./tests_tmp/1.bit_matrix");
    vmw.load(matl);
    EXPECT_TRUE(!memcmp(mat.matrix, rev->matrix, 32));
  }
  {
    VectorMatrixWriter vmw("./tests_tmp/1.bit_matrix.lz4", 0, 0, 0, 0, 1, false);
    vmw.dump(mat);
  }
  {
    BitMatrix matl(16, 2, true);
    VectorMatrixReader vmw("./tests_tmp/1.bit_matrix.lz4");
    vmw.load(matl);
    EXPECT_TRUE(!memcmp(mat.matrix, rev->matrix, 32));
  }
  delete trp;
  delete rev;
}
