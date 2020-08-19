#include "kmtricks/sequences.hpp"

using namespace km;

// set a custom hasher
template<typename K>
class CustomHasher : public Hasher<K>
{
public:
  CustomHasher() {};
  ~CustomHasher() {};
  uint64_t operator() (K kmer, uint64_t seed) override
  {
    return 1;
  }
};

// set a custom validator for minimizers
template<typename K>
class CustomValidator : public Validator<K>
{
public:
  CustomValidator() {};
  ~CustomValidator() {};
  bool operator() (K value, size_t size) override
  {
    if (value > 0) return true;
    return false;
  }
};

int main(int argc, char* argv[])
{
  // Kmer common
  Kmer<uint64_t> kmer1("ACGCTCTTTTT", false); // kmer, canonical
  cout << "k-mer size " << to_string(kmer1.size()) << endl;
  cout << "k-mer value " << to_string(kmer1.value()) << endl;
  cout << "k-mer str value " << kmer1.str_value() << endl;
  cout << "k-mer rev value " << to_string(kmer1.rev_comp()) << endl;
  cout << "k-mer str rev " << kmer1.str_rev_comp() << endl;


  // Kmer hash
  cout << "k-mer hash " << to_string(kmer1.hash()) << endl;


  CustomHasher<uint64_t> *custom = new CustomHasher<uint64_t>(); // custom hasher from Hasher<K>
  kmer1.set_hasher(custom); // set a custom hasher
  cout << "k-mer hash custom " << to_string(kmer1.hash()) << endl; // or use test.hash(custom)
  delete custom;

  // change to canonical mode
  kmer1.use_canonical();
  cout << "use_canonical() -> k-mer str value " << kmer1.str_value() << endl;

  // Canonical  at construct time
  Kmer<uint64_t> kmer2("ACGCTCTTTTT", true);
  cout << "Canonical form at construct: str value " << kmer2.str_value() << endl;

  //using specific encoding scheme
  cout << endl << "Custom encoding" << endl;
  cout << "default encoding, k-mer value " << to_string(kmer2.value()) << ", str value " << kmer2.str_value() << endl;

  uchar mymap[4] = { 'T', 'A', 'C', 'G' };
  Code<uint64_t> my_encoder(mymap);
  Kmer<uint64_t> kmer3("ACGCTCTTTTT", true, &my_encoder);
  cout << "With custom encoding, k-mer value " << to_string(kmer3.value()) <<", str value " <<  kmer3.str_value() << endl;

  // Superk

  Superk<uint64_t> superk("AGCAGAGCAAAAGAAAAAGAAAACGAGAAAAACAAAGACAACGAAACTTATAATTTATATCACTACGATTATAAAAAAACTTATTATATTTAAT", 31);
  Kmer<uint64_t> firstk = superk.get_first();
  cout << "first : " << firstk.str_value() << endl;
  cout << "Superk: " << superk.str_value() << endl;
  cout << "Size: " << to_string(superk.size()) << endl;

  for (size_t i=0; i<superk.nb_kmers(); i++)
    cout << to_string(i) << " " << superk.get_kmer(i, false).str_value() << endl;

  // get kmer at pos N
  Kmer<uint64_t> kmerN = superk.get_kmer(0, false);
  cout << "kmer pos 0: " << kmerN.str_value() << endl;

  // or use :
  Kmer<uint64_t> kmerX(false); // true for canonical form
  // if you use custom encoding, you can use encoding from superk:
  // Kmer<uint64_t> kmerCustom(false, superk.get_encoding());
  superk.get_kmer(1, &kmerX);
  cout << "Kmer pos 1: " << kmerX.str_value() << endl;
  superk.get_kmer(2, &kmerX);
  cout << "Kmer pos 2: " << kmerX.str_value() << endl;
  superk.get_kmer(3, &kmerX);
  cout << "Kmer pos 3: " << kmerX.str_value() << endl;
  superk.get_kmer(12, &kmerX);
  cout << "Kmer pos 12: " << kmerX.str_value() << endl;


  // Minimizer
  // From kmer
  cout << endl << "Minim from kmer" << endl;
  Kmer<uint64_t> kmer4("ACGCTCTTTTT", true);
  Minimizer<uint64_t> miniKmer(&kmer4, 10, true); // build minimizer of size 10, with check_validity = true
  cout << "miniKmer str: " << miniKmer.str_value() << " miniKmer val: " << to_string(miniKmer.value()) << endl;
  cout << "Validator returns false, so returns default minimizer, without validation :" << endl;
  Minimizer<uint64_t> miniKmerNv(&kmer4, 10, false);
  cout << "miniKmer str: " << miniKmerNv.str_value() << " miniKmer val: " << to_string(miniKmerNv.value()) << endl;
  cout << "Use custom default minim" << endl;
  Minimizer<uint64_t> miniKmerD(&kmer4, 10, true, 500);
  cout << "miniKmer str: " << miniKmerD.str_value() << " miniKmer val: " << to_string(miniKmerD.value()) << endl;
  //From superk
  cout << endl << "Minim from superk" << endl;
  Minimizer<uint64_t> miniSuperk(&superk, 10, true);
  cout << "miniSuperk str: " << miniSuperk.str_value() << " miniSuperk val: " << to_string(miniSuperk.value()) << endl;

  // use custom validator
  CustomValidator<uint64_t> *customValid = new CustomValidator<uint64_t>();
  Minimizer<uint64_t> miniKmer2(&kmer4, 10, customValid);
  // custom valid and custom default
  uint64_t def_minim = 500;
  Minimizer<uint64_t> miniKmer3(&kmer4, 10, customValid, def_minim);

  delete customValid;
  return 0;
}