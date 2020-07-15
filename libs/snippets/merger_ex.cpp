#include "kmtricks/merger.hpp"
#include <iostream>

int main(int argc, char* argv[])
{
  try
  {
    typedef selectK<31>::type   KType;           // Size of k-mers
    typedef selectC<255>::type  CType;           // Max count

    string fof = argv[1];
    Merger<KType, CType> m(fof, 1, 1, 0, true); // fof, min abundance, min recurrence, header size (0 for headerless file)
    while(!m.end)
    {
      m.next();
      if ( m.keep )
      {
        cout << to_string(m.m_khash);           // ascii int value (hash or kmer)
        //cout << m.getk(31, bToN);             // ascii kmer
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
