#include <kmtricks/plugin.hpp>

// DMAX_C is a compile definition set by cmake
using count_type = typename km::selectC<DMAX_C>::type;


class BasicEx : public km::IMergePlugin
{
public:
  BasicEx() = default;
private:
  unsigned int m_threshold {0};

public:
  // Override process_kmer
  // Discard lines which contain abundances less than a threshold
  bool process_kmer(const uint64_t* kmer_data, std::vector<count_type>& count_vector) override
  {
    for (auto& c : count_vector)
      if (c < m_threshold)
        return false;
    return true;
  }

  // Override configure (not necessary if you don't need configuration)
  // The string is passed to kmtricks with --plugin-config
  // Here it's a simple example where the string is a threshold
  // It could be a path to a config file for instance
  void configure(const std::string& s) override
  {
    m_threshold = std::stoll(s);
  }
};

// Make the plugin loadable
extern "C" std::string plugin_name() { return "BasicEx"; }
extern "C" int use_template() { return 0; }
extern "C" km::IMergePlugin* create0() { return new BasicEx(); }
extern "C" void destroy(km::IMergePlugin* p) { delete p; }

