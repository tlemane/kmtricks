
#include <kmtricks/plugin.hpp>


using count_type = typename km::selectC<DMAX_C>::type;

template<size_t MAX_K>
class ExamplePlugin : public km::IMergePlugin
{
public:
  ExamplePlugin() = default;

  bool process_kmer(const uint64_t* kmer_data, std::vector<count_type>& count_vector) override
  {
    m_kmer.set64_p(kmer_data);

    if (m_kmer.at(0) == 'A') // set count to 42 for k-mers start with A
    {
      for (auto& c : count_vector)
        c = 42;
    }
    if (m_kmer.at(0) == 'C') // discard line when k-mers start with C
      return false;
    return true;
  }

  bool process_hash(uint64_t h, std::vector<count_type>& count_vector) override
  {
    return true;
  }

  void set_kmer_size(size_t kmer_size) override
  {
    m_kmer.set_k(kmer_size);
  }

private:
  km::Kmer<MAX_K> m_kmer;
};


extern "C" std::string plugin_name() { return "ExamplePluginTemplate"; }
extern "C" int use_template() { return 1; }
extern "C" km::IMergePlugin* create32() { return new ExamplePlugin<32>; }
extern "C" km::IMergePlugin* create64() { return new ExamplePlugin<64>; }

// enables kmtricks usage with max --kmer-size 63
// Depending on --kmer-size, kmtricks will check if the right createX() function is available
// before start computation.

extern "C" void destroy(km::IMergePlugin* p) { delete p; }
