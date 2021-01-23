#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <cassert>
#include <unordered_map>
#include <kmtricks/io.hpp>
#include <kmtricks/sequences.hpp>
using namespace km;
using namespace std;

random_device rd;
mt19937 mt(rd());
uniform_int_distribution<uint64_t> distK;
uniform_int_distribution<uint16_t> distC;

typedef uint64_t kt;
typedef uint8_t ct;

vector<uint16_t> get_rand_count(int n)
{
  vector<uint16_t> v(n);
  for (int i=0; i<n; i++)
    v[i] = distC(mt)%64;
  return v;
}

int main(int argc, char* argv[])
{
  Kmer<kt>   kmer("ACGT"); // a k-mer
  ct         count(10); // a count
  vector<ct> counts_vec = {1, 2, 3, 4, 5, 6}; // a count vector
  vector<char> B(1); // bit vector, e.g. p/a vector
  auto nbbits = 8;
  Kmer<kt>   inkmer(false);
  ct         incount;
  vector<ct> incounts;
  vector<char> inb(1);
  auto fileid = 0;
  auto partid = 0;
  auto kmer_size = 4;
  auto is_hash = 0;
  auto is_compressed = 0;
  
  {
    KmerFile<OUT, kt, ct> kmer_file(
      "KmerFile.kmer", fileid, partid, kmer_size, is_hash, is_compressed);
    kmer_file.write(kmer, count);
  }
  {
    KmerFile<IN, kt, ct> kmer_file(
      "KmerFile.kmer");
    kmer_file.read(inkmer, incount);
  }

  {
    CountMatrixFile<OUT, kt, ct, matrix_t::BIN> cmat(
      "CountMatrix.mat", partid, 6, kmer_size, is_hash, is_compressed);
    cmat.write(kmer, counts_vec);
  }
  {
    CountMatrixFile<IN, kt, ct, matrix_t::BIN> cmat(
      "CountMatrix.mat");
    cmat.read(inkmer, incounts);
  }

  {
    PAMatrixFile<OUT, kt> pam(
      "PaMat.pa", partid, 6, kmer_size, is_hash, is_compressed);
    pam.write(kmer, B);
  }
  {
    PAMatrixFile<IN, kt> pam(
      "PaMat.pa");
    pam.read(inkmer, inb);
  }

  {
    BitVectorFile<OUT> bv(
      "BitVec.vec", fileid, partid, nbbits, is_compressed);
    bv.write(B);
  }
  {
    BitVectorFile<IN> bv(
      "BitVec.vec");
    bv.write(inb);
  }

  {
    BitMatrixFile<OUT, matrix_t::BF> bmat(
      "BitMat.mat", partid, 8, 8, is_compressed);
    bmat.write(B);
  }
  {
    BitMatrixFile<IN, matrix_t::BF> bmat(
      "BitMat.mat");
    bmat.read(inb);
  }
  return 0;
}