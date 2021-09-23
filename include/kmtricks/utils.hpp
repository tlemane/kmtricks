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

#pragma once
#include <string>
#include <fstream>
#include <random>
#include <sstream>
#include <climits>
#include <thread>
#include <condition_variable>
#include <queue>
#include <vector>
#include <unistd.h>
#include <functional>
#include <filesystem>
#include <cmath>
#include <cstdlib>
#include <sys/resource.h>

namespace fs = std::filesystem;

#if __APPLE__
#include <mach-o/dyld.h>
#include <mach/mach.h>
#endif

#ifndef KMTRICKS_PUBLIC
#include <spdlog/spdlog.h>
#include <fmt/format.h>
#endif

#include <kmtricks/exceptions.hpp>

#define ROUND_UP(n, m) ((n + m - 1) / m) * m
#define NBYTES(m) (m + 7) / 8
#define BITMASK(b) (1 << ((b) % CHAR_BIT))
#define BITSLOT(b) ((b) / CHAR_BIT)
#define BITSET(a, b) ((a)[BITSLOT(b)] |= BITMASK(b))
#define BITCHECK(a, b) ((a)[BITSLOT(b)] & BITMASK(b))

#if defined(__GNUC__) || defined(__clang__)
#define KM_LIKELY(x) __builtin_expect((x), 1)
#define KM_UNLIKELY(x) __builtin_expect((x), 0)
#else
#define KM_LIKELY(x) x
#define KM_UNLIKELY(x) x
#endif

namespace km {

#ifndef KMTRICKS_PUBLIC
enum class VerbosityLevel
{
  DEBUG,
  INFO,
  WARNING,
  ERROR
};

void set_verbosity_level(const std::string& level);
VerbosityLevel str_to_verbosity_level(const std::string& str_level);

template<typename T>
void check_fstream_good(const std::string& path, const T& stream)
{
  if (!stream.good())
  {
    if constexpr(std::is_same_v<T, std::ofstream>)
      throw IOError(fmt::format("Unable to write at {}.", path));
    else
      throw IOError(fmt::format("Unable to read at {}.", path));
  }
}
#else
template<typename T>
void check_fstream_good(const std::string& path, const T& stream)
{
  if (!stream.good())
  {
    if constexpr(std::is_same_v<T, std::ofstream>)
      throw IOError("Unable to write at " + path);
    else
      throw IOError("Unable to read at " + path);
  }
}
#endif

template<typename T>
void set_bit_vector(std::vector<uint8_t>& bit_vec, const std::vector<T>& count_vec)
{
  assert(count_vec.size() <= bit_vec.size()*8);
  std::fill(bit_vec.begin(), bit_vec.end(), 0);
  for (size_t i=0; i<count_vec.size(); i++)
  {
    if (count_vec[i])
    {
      BITSET(bit_vec, i);
    }
  }
}

template<size_t MAX_K>
uint64_t get_required_memory(size_t nb_kmers)
{
  return nb_kmers * (((MAX_K + 31) / 32) * 8) + 8192;
}

template<size_t MAX_K>
uint64_t get_required_memory_hash(size_t nb_kmers)
{
  return nb_kmers * (sizeof(uint64_t)) + 8192;
}

inline std::string get_uname_sr()
{
  std::array<char, 256> buffer;
  std::string result;
  auto pipe = popen("uname -sr", "r");

  if (!pipe)
    throw std::runtime_error("popen failed.");

  while (!feof(pipe))
  {
    if (fgets(buffer.data(), 256, pipe) != nullptr) result += buffer.data();
  }
  auto rc = pclose(pipe);
  if (rc != 0)
    throw std::runtime_error("pclose failed.");

  return result;
}

inline std::string random_dna_seq(size_t size)
{
  static const char alpha[] = "ACGT";
  static std::default_random_engine g;
  static std::uniform_int_distribution<int> dist(0, 3);
  std::string seq;
  for (size_t i = 0; i < size; i++)
    seq.push_back(alpha[dist(g)]);
  return seq;
}

template<typename T>
inline std::vector<T> random_count_vector(size_t size)
{
  static std::default_random_engine g;
  static std::uniform_int_distribution<int> dist(0, std::numeric_limits<T>::max());
  std::vector<T> vec;
  for (size_t i = 0; i < size; i++)
    vec.push_back(dist(g));
  return vec;
}

template<typename T>
inline std::string vec_to_str(const std::vector<T>& vec)
{
  std::stringstream ss;
  for (const auto& v : vec)
    ss << std::to_string(v) << " ";
  return ss.str();
}

inline std::vector<char*> str_vec_to_c_str(const std::vector<std::string>& vec)
{
  std::vector<char*> cstrings(vec.size());
  for (auto& s : vec)
    cstrings.push_back(const_cast<char*>(s.c_str()));
  return cstrings;
}

template<typename T>
inline std::ostream& write_vector(std::ostream& out, std::vector<T> vec, char delim = ' ')
{
  for (auto& e: vec)
    out << std::to_string(e) << delim;
  return out;
}

inline size_t get_peak_rss()
{
  struct rusage rusage;
  getrusage(RUSAGE_SELF, &rusage);
#if __APPLE__
  return static_cast<size_t>(rusage.ru_maxrss / 1024L);
#else
  return static_cast<size_t>(rusage.ru_maxrss);
#endif
}

inline size_t get_current_rss()
{
#if __APPLE__
  struct mach_task_basic_info info;
  mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
  if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
    (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
    return (size_t)0L;
  return (size_t)info.resident_size;
#else
  long rss = 0;
  FILE *fp = NULL;
  if ((fp = fopen("/proc/self/statm", "r")) == NULL)
    return 0;
  if (fscanf(fp, "%*s%ld", &rss) != 1)
  {
    fclose(fp);
    return 0;
  }
  fclose(fp);
  return rss * sysconf(_SC_PAGE_SIZE);
#endif
}

inline std::tuple<int64_t, int64_t> get_prlimit_nofile()
{
  struct rlimit rlim;
  getrlimit(RLIMIT_NOFILE, &rlim);
  return std::make_tuple(rlim.rlim_cur, rlim.rlim_max);
}

inline double bloom_fp(size_t m, size_t n, size_t k = 1)
{
  static double e = std::exp(1.0);
  return std::pow(1.0 - std::pow(e, (-(k * static_cast<double>(n)) / static_cast<double>(m))), static_cast<double>(k));
}

inline double bloom_estimate(size_t m, size_t k, size_t x)
{
  return -(static_cast<double>(m)/static_cast<double>(k) * std::log(1.0 - static_cast<double>(x)/static_cast<double>(m)));

}
class Eraser
{
public:
  Eraser(size_t nb_threads = 1)
  {
    for (size_t i=0; i<nb_threads; i++)
      m_pool.push_back(std::thread(&Eraser::worker, this));
  }

  static Eraser& get()
  {
    static Eraser singleton;
    return singleton;
  }

  void erase(const std::string& path)
  {
    {
      std::unique_lock<std::mutex> lock(m_mutex);
      m_queue.push(path);
    }
    m_condition.notify_one();
  }

  void join()
  {
    {
      std::unique_lock<std::mutex> lock(m_mutex);
      m_stop = true;
    }
    m_condition.notify_all();
    for (auto& t : m_pool)
      if (t.joinable())
        t.join();
  }
private:
  void worker()
  {
    while (true)
    {
      std::string path;
      {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_condition.wait(lock, [this]{return this->m_stop || !this->m_queue.empty();});
        if (m_stop && m_queue.empty())
          return;
        path = std::move(m_queue.front());
        m_queue.pop();
      }
      fs::remove(path);
    }
  }

private:
  std::vector<std::thread> m_pool;
  std::mutex m_mutex;
  std::condition_variable m_condition;
  std::queue<std::string> m_queue;
  bool m_stop {false};
};

template<size_t C>
struct requiredC
{
  enum
  {
    value = C <= 0xFF ? 8 : C <= 0xFFFF ? 16 : 32
  };
};

template<int bits> struct select_;

template<> struct select_ <8> {typedef uint8_t type;};
template<> struct select_ <16> {typedef uint16_t type;};
template<> struct select_ <32> {typedef uint32_t type;};

template<size_t C>
struct selectC : select_<requiredC<C>::value> {};

};
