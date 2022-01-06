#include <cmath>
#include <map>
#include <kff_io.hpp>
#include <string>
#include <vector>
#include <tuple>
#include <memory>

#include "encoding.hpp"

using uint = unsigned int;

using ptr_t = std::unique_ptr<std::uint8_t[]>;
using bucket_t = std::tuple<ptr_t, ptr_t, uint, uint>;

int main(int argc, char* argv[])
{
  std::vector<std::string> filenames(argv + 1, argv + argc);

  Kff_file infile(filenames[0], "r");
  Section_GV gv(&infile);
  uint msize = infile.global_vars["m"];
  uint kmer_size = infile.global_vars["k"];
  gv.close();
  infile.close();

  Stringifyer st(infile.encoding);

  std::unordered_map<std::string, std::vector<bucket_t>> data; //(std::pow(4, msize));

  std::size_t counter {0};
  for (auto& file : filenames)
  {
    Kff_file infile(file, "r");

    char section_type = infile.read_section_type();

    while (infile.tellp() < infile.end_position)
    {

      if (section_type == 'i')
      {
        Section_Index si(&infile);
        si.close();
      }
      if (section_type == 'v')
      {
        Section_GV sgv(&infile);
        sgv.close();
      }
      else if (section_type == 'm')
      {
        Section_Minimizer sm(&infile);
        uint k = infile.global_vars["k"];
        uint m = infile.global_vars["m"];
        uint max = infile.global_vars["max"];
        uint data_size = infile.global_vars["data_size"];

        counter++;

        uint max_nucl = k - m + max - 1;


        for (uint i=0; i<sm.nb_blocks; i++)
        {
          std::unique_ptr<std::uint8_t[]> seq(new std::uint8_t[max_nucl / + 1]);
          std::unique_ptr<std::uint8_t[]> data_bytes(new std::uint8_t[data_size * max]);
          std::uint64_t mini_pos;
          uint nb_kmers = sm.read_compacted_sequence_without_mini(seq.get(), data_bytes.get(), mini_pos);
          data[st.translate(sm.minimizer, m)].push_back(std::make_tuple(std::move(seq), std::move(data_bytes), mini_pos, nb_kmers));
        }

      }
      section_type = infile.read_section_type();
    }
    infile.close();
  }

  Kff_file outfile("out.kff", "w");
  std::uint8_t encoding[] = {0, 1, 3, 2};
  outfile.write_encoding(encoding);

  Section_GV sgv(&outfile);
  sgv.write_var("k", kmer_size);
  sgv.write_var("m", msize);
  sgv.write_var("max", kmer_size - msize + 1);
  sgv.write_var("data_size", 1);
  sgv.close();

  Binarizer bin(encoding);

  for (auto& [m, v] : data)
  {
    Section_Minimizer sm(&outfile);
    std::uint8_t minim[m.size()];
    bin.translate(m, msize, minim);
    sm.write_minimizer(minim);

    for (auto& [i, j, k, l] : v)
    {
      sm.write_compacted_sequence_without_mini(i.get(), kmer_size-msize+l, k, j.get());
    }

    sm.close();
  }

  outfile.close();

  std::cerr << "Before merge : " << counter  << ", After merge :  " << data.size() << std::endl;
}
