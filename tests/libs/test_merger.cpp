#include <criterion/criterion.h>
#include <kmtricks/merger.hpp>

using namespace km;
typedef uint64_t kt;
typedef uint8_t  kc;

Test(merger, merger_build)
{
  string path = "./merge_data/fof.txt";
  Merger<kt, kc> m(path, 1, 1, 0, true);
  cr_assert(m.nb_files == 2);
}

Test(merger, merger_next)
{
  string path = "./merge_data/fof.txt";
  Merger<kt, kc> m(path, 1, 2, 0, true);

  kt values[5] = {0, 134, 234, 300, 302};
  kc count1[5] = {1, 31, 1, 8, 12};
  kc count2[5] = {3, 6, 100, 2, 1};

  int i = 0;
  while(!m.end)
  {
    m.next();
    if (m.keep)
    {
      cr_assert(m.m_khash == values[i]);
      cr_assert(m.counts[0] == count1[i]);
      cr_assert(m.counts[1] == count2[i]);
      i++;
    }
  }
}

Test(merger, merger_get_kmer)
{
  string path = "./merge_data/fof.txt";
  Merger<kt, kc> m(path, 1, 2, 0, true);

  kt values[5] = {0, 134, 234, 300, 302};

  int i = 0;
  while(!m.end)
  {
    m.next();
    if (m.keep)
    {
      Kmer<kt> k(values[i], 32, false);
      Kmer<kt> k2 = m.get_kmer(32);
      cr_assert(k == k2);
      i++;
    }
  }
}