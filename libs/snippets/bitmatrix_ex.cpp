/*****************************************************************************
 *   kmtricks
 *   Authors: T. Lemane
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#include "kmtricks/bitmatrix.hpp"
using namespace km;
int main(int argc, char* argv[])
{
  try
  {
    int nb = 16; // nb line in bits
    int m = 2;   // nb cols in bytes
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
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
  }
  return 0;
}
