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

#include "kmtricks/minim_repart.hpp"

#include <string>
#include <iostream>

int main(int argc, char* argv[])
{
  string path = argv[1];
  RepartFile f(path);
  string kmer = "GAGCAGCACAAACGAGACAC";
  // string kmer = "AAAAAAAAAAAAAAAAAAAAAAA"; // no valid minimizer
  size_t sk = kmer.size();
  size_t sm = 10;
  MinimRepart<uint64_t> m(f);

  uint64_t kmerbin = m.seq_to_int(kmer, sk);
  uint64_t revcomp = m.rev_comp(kmerbin, sk);
  string revcomp_str = m.int_to_str(revcomp, sk);



  cout << "KMER : " <<  kmer << " " << to_string(kmerbin) << endl;
  cout << "REV  : " << revcomp_str << " " << to_string(revcomp) << endl;

  uint64_t minim = m.get_minim_from_str(kmer, sk, sm);
  if (minim == numeric_limits<uint64_t>::max()){
    cerr << "No minimizer found" << endl;
    return 1;
  }
  string   minim_str = m.int_to_str(minim, sm);

  cout << "MINIM : " << minim_str << " " << to_string(minim) << endl;

  uint part = m.get_partition(minim);

  cout << "Part = " << part << endl;
}