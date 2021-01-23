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

#include "kmtricks/merger.hpp"
#include <iostream>

using namespace km;

int main(int argc, char* argv[])
{
  try
  {
    typedef selectK<31>::type   KType;           // Size of k-mers
    typedef selectC<255>::type  CType;           // Max count

    string fof = argv[1];
    int amin = atoi(argv[2]);
    int rmin = atoi(argv[3]);
    Merger<KType, CType, KmerFile<IN, KType, CType>> m(fof, amin, rmin, 0, true); // fof, min abundance, min recurrence, header size (0 for headerless file)
    while(!m.end)
    {
      m.next();
      if ( m.keep )
      {
        cout << to_string(m.m_khash);           // ascii int value (hash or kmer)
        //cout << m.get_kmer(31).str_value();   // get_kmer(size) return Kmer<KType>
        for (size_t i=0; i<m.nb_files; i++)     // counts
        {
          cout << " " << to_string(m.counts[i]);
        }
        cout << "\n";
      }
    }
  }
  catch ( exception& e)
  {
    cout << e.what() << endl;
  }
  return 0;
}
