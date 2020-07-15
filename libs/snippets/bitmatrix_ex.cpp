#include "kmtricks/bitmatrix.hpp"

int main(int argc, char* argv[])
{
  int nb = 16; // nb line in bits
  int m = 1;   // nb cols in bytes
  BitMatrix* mat = new BitMatrix(nb, m, true);

  for (int i=0; i<nb; i++)
    for (int j=0; j<m; j++)
    {
      mat->set_bit(i, 0, 1);
    }
  mat->set_bit(1,0,0);
  BitMatrix* trp = mat->transpose();

  mat->print_bytes();
  trp->print_bytes();

  mat->print_bits();
  trp->print_bits();

  trp->dump("test_trp.mat");
}
