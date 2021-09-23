#include <kmtricks/plugin.hpp>
#include <spdlog/spdlog.h>

using count_type = typename km::selectC<DMAX_C>::type;

class ExamplePlugin : public km::IMergePlugin
{
public:
  ExamplePlugin() = default;

  bool process_kmer(const uint64_t* kmer_data, std::vector<count_type>& count_vector) override
  {
    for (auto& c : count_vector) // set all 0 to 42
      if (c == 0)
        c = 42;
    return true;
  }

  bool process_hash(uint64_t h, std::vector<count_type>& count_vector) override
  {
    return true;
  }
};


extern "C" std::string plugin_name() { return "ExamplePlugin"; }
extern "C" int use_template() { return 0; }
extern "C" km::IMergePlugin* create0() { return new ExamplePlugin(); }
extern "C" void destroy(km::IMergePlugin* p) { delete p; }
