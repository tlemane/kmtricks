#include <criterion/criterion.h>
#include "skreader.hpp"
#include "sequences.hpp"

using namespace km;
typedef uint64_t kt;

Test(skstorage, skstorage_build)
{
  string path = "./skreader_data/sk_part";
  string prefix = "superKparts.";
  SuperkStorage s(path, prefix, 4);

  cr_assert(s.nb_files() == 4);

  for (int i=0; i<s.nb_files(); i++)
  {
    cr_assert(s._parts[i].is_open());
  }
}

Test(skstorage, skstorage_read_block)
{
  string path = "./skreader_data/sk_part";
  string prefix = "superKparts.";
  SuperkStorage s(path, prefix, 4);

  uchar *buffer = 0;
  uint  size = 0;
  uint  nb_bytes = 0;

  int r = s.read_block(&buffer, &size, &nb_bytes, 0);
  cr_assert(r == 72);
  cr_assert(*buffer == 4);

  Superk<kt> superk(++buffer, 23, 20, true);
  cr_assert(superk.str_value() == "CATACAGAGACAGCAGCAGAGCA");
}

Test(skreader, skreader_build)
{
  string path = "./skreader_data/sk_part";
  string prefix = "superKparts.";
  SuperkStorage s(path, prefix, 4);
  size_t kmersize = 20;
  SuperkReader<kt> reader(&s, kmersize);
  Superk<kt> superk(kmersize);
}

Test(skreader, skreader_next)
{
  string path = "./skreader_data/sk_part";
  string prefix = "superKparts.";
  SuperkStorage s(path, prefix, 4);
  size_t kmersize = 20;
  SuperkReader<kt> reader(&s, kmersize);
  Superk<kt> superk(kmersize);
  reader.next_superk(0, &superk);
  cr_assert(superk.str_value() == "CATACAGAGACAGCAGCAGAGCA");
  reader.next_superk(0, &superk);
  cr_assert(superk.str_value() == "GAGCAGCACAAACGAGACACAAAAAAAGAG");
}