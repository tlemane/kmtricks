#include <criterion/criterion.h>
#include <time.h>
#include <kmtricks/bitmatrix.hpp>

typedef uint64_t kt;
static uchar zero[32];
static uchar one[32];

using namespace km;

void setup(void)
{
  memset(zero, 0, 32);
  memset(one, 0xFF, 32);
}

Test(BitMatrix, bitmatrix_build, .init = setup)
{
  BitMatrix mat(16, 2, true, false); // 16 line in bits, 2 cols in bytes -> 16*16
  cr_assert(!memcmp(mat.matrix, zero, 32));
  BitMatrix mat1(16, 2, true, true);
  cr_assert(!memcmp(mat1.matrix, one, 32));
}

Test(BitMatrix, bitmatrix_clear, .init = setup)
{
 BitMatrix mat(16, 2, true, true);
 mat.clear();
 cr_assert(!(memcmp(mat.matrix, zero, 32)));
}

Test(BitMatrix, bitmatrix_set_bit)
{
  BitMatrix mat(16, 2, true); // 16 line in bits, 2 cols in bytes -> 16*16
  cr_assert(mat.get_bit(4, 6) == 0);
  mat.set_bit(4, 6, true);
  cr_assert(mat.get_bit(4, 6) == 1);
}

Test(BitMatrix, bitmatrix_tog_bit)
{
  BitMatrix mat(16, 2, true); // 16 line in bits, 2 cols in bytes -> 16*16
  mat.set_bit(4, 6, true);
  cr_assert(mat.get_bit(4, 6) == 1);
  mat.tog_bit(4, 6);
  cr_assert(mat.get_bit(4, 6) == 0);
}

Test(BitMatrix, bitmatrix_set_byte)
{
  BitMatrix mat(16, 2, true);
  mat.set_byte(1, 1, 0x80);
  cr_assert(mat.get_bit(1,15) == true);
}

Test(BitMatrix, bitmatrix_tog_byte)
{
  BitMatrix mat(16, 2, true);
  mat.set_byte(1, 1, 0x80);
  mat.tog_byte(1, 1);
  cr_assert(mat.get_byte(1, 1) == 0x7F);
  cr_assert(mat.get_bit(1,15) == false);
}

Test(BitMatrix, bitmatrix_tranpose)
{
  srand(time(NULL));
  int i, j;
  BitMatrix mat(16, 2, true);
  for (int n=0; n<20; n++)
  {
    i = rand() % 16;
    j = rand() % 16;
    mat.set_bit(i, j, true);
  }

  BitMatrix *trp = mat.transpose();
  BitMatrix *rev = trp->transpose();

  cr_assert(!memcmp(mat.matrix, rev->matrix, 32));

  delete trp;
  delete rev;
}
