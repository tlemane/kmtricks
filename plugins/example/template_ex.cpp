#include <kmtricks/plugin.hpp>

// Same as BasicEx
using count_type = typename km::selectC<DMAX_C>::type;

template<std::size_t MAX_K>
class TemplateEx : public km::IMergePlugin
{
public:
  TemplateEx() = default;
private:
  unsigned int m_threshold {0};
  // Declare a k-mer
  km::Kmer<MAX_K> m_kmer;

public:
  // same as BasicEx
  void configure(const std::string& s)
  {
    m_threshold = std::stoll(s);
  }

  // Override set_kmer_size to pass the k-mer size to m_kmer
  void set_kmer_size(size_t kmer_size) override
  {
    this->m_kmer_size = kmer_size;
    m_kmer.set_k(this->m_kmer_size);
  }

  // Override process_kmer
  // Discard lines which contain abundances less than a threshold if the k-mer starts with 'A'
  bool process_kmer(const uint64_t* kmer_data, std::vector<count_type>& count_vector) override
  {
    m_kmer.set64_p(kmer_data);
    if (m_kmer.at(0) == 'A')
    {
      for (auto& c : count_vector)
      {
        if (c < m_threshold)
        {
          return false;
        }
      }
    }
    return true;
  }
};

// Make the plugin loadable
extern "C" std::string plugin_name() { return "TemplateEx"; }
extern "C" int use_template() { return 1; }
extern "C" km::IMergePlugin* create32() { return new TemplateEx<32>(); } // call if --kmer-size < 32
extern "C" km::IMergePlugin* create64() { return new TemplateEx<64>(); } // call if --kmer-size < 64
extern "C" km::IMergePlugin* create512() { return new TemplateEx<512>(); } // call if --kmer-size < 512

// With create32, create64 and create512, the plugin supports k-mer size in [8, 64) and [480, 512)j
