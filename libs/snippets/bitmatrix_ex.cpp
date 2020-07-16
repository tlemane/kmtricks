#include "kmtricks/bitmatrix.hpp"

int main(int argc, char* argv[])
{
  int nb = 32; // nb line in bits
  int m = 8;   // nb cols in bytes
  BitMatrix* mat = new BitMatrix(nb, m, true);

  for (int i=0; i<nb; i++)
    for (int j=0; j<m; j++)
    {
      mat->set_bit(i, 0, 1);
    }
  mat->set_bit(1,0,0);
  BitMatrix* trp = mat->transpose();
  BitMatrix* trp2 = trp->transpose();

  for (int i=0; i<nb*m; i++)
    assert(mat->matrix[i] == trp2->matrix[i]);

  mat->print_bits();
  trp->print_bits();

  delete mat;
  delete trp;
  delete trp2;

  return 0;
}
