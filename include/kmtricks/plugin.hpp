#pragma once
#include <string>
#include <vector>
#include <cstdint>
#include <kmtricks/kmer.hpp>

#define KMTRICKS_PUBLIC
#include <kmtricks/utils.hpp>

namespace km {

class IMergePlugin
{
public:
  IMergePlugin() = default;
  virtual ~IMergePlugin() {}
  virtual void set_out_dir(const std::string& s) final { m_output_directory = s; }
  virtual void set_partition(size_t p) final { m_partition = p; }
  virtual void set_kmer_size(const size_t kmer_size) { m_kmer_size = kmer_size; }

  virtual void configure(const std::string& s) {}

  virtual bool process_kmer(const uint64_t* kmer_data, std::vector<typename selectC<DMAX_C>::type>& count_vector) { return true; }
  virtual bool process_hash(uint64_t h, std::vector<typename selectC<DMAX_C>::type>& count_vector) { return true; }

protected:
  std::string m_output_directory;
  size_t m_kmer_size;
  size_t m_partition;
};

}
