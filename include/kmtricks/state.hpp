#ifndef STATE_HPP_OSPLYBCL
#define STATE_HPP_OSPLYBCL

#include <spdlog/spdlog.h>
#include <mutex>
#include <vector>
#include <fstream>
#include <kmtricks/spinlock.hpp>
#include <cstdint>

namespace km {

  class state
  {
    public:
      static state& get()
      {
        static state s;
        return s;
      }

    public:

      void init(const std::string path, std::size_t f, std::size_t p)
      {
        m_path = path;
        m_nb = f;
        m_part = p;
        m_superk.resize(f, 0);
        m_count.resize(f * p, 0);
        m_merge.resize(p, 0);
      }

      void init_f(const std::string& f, std::size_t n, std::size_t p)
      {
        m_nb = n;
        m_part = p;
        std::ifstream inf(f, std::ios::in | std::ios::binary);

        m_superk.resize(m_nb, 0);
        m_count.resize(m_nb * m_part, 0);
        m_merge.resize(m_part, 0);

        inf.read(reinterpret_cast<char*>(&m_repart), sizeof(m_repart));
        inf.read(reinterpret_cast<char*>(&m_config), sizeof(m_config));
        inf.read(reinterpret_cast<char*>(m_superk.data()), sizeof(std::uint8_t) * n);
        inf.read(reinterpret_cast<char*>(m_count.data()), sizeof(std::uint8_t) * (n*p));
        inf.read(reinterpret_cast<char*>(m_merge.data()), sizeof(std::uint8_t) * p);
      }

      void write()
      {
        std::unique_lock<spinlock> _1(m_lock_1);
        std::unique_lock<spinlock> _2(m_lock_2);
        std::ofstream inf(m_path, std::ios::out | std::ios::binary);

        inf.write(reinterpret_cast<char*>(&m_config), sizeof(m_config));
        inf.write(reinterpret_cast<char*>(&m_repart), sizeof(m_repart));
        inf.write(reinterpret_cast<char*>(m_superk.data()), sizeof(std::uint8_t) * m_nb);
        inf.write(reinterpret_cast<char*>(m_count.data()), sizeof(std::uint8_t) * (m_nb*m_part));
        inf.write(reinterpret_cast<char*>(m_merge.data()), sizeof(std::uint8_t) * m_part);
      }

    public:
      bool repart() { return m_repart; }
      bool config() { return m_config; }
      bool superk(std::size_t f) { return m_superk[f]; }
      bool count(std::size_t f, std::size_t p) {
        return m_count[p + m_part * f];
      }
        bool merge(std::size_t p) { return m_merge[p]; }

      void repart_done() { m_repart = true; }
      void config_done() { m_config = true; }

      void superk_done(std::size_t f)
      {
        std::unique_lock<spinlock> _(m_lock_1);
        m_superk[f] = 1;
      }

      void count_done(std::size_t f, std::size_t p)
      {
        std::unique_lock<spinlock> _(m_lock_2);
        m_count[p + m_part * f] = 1;
      }

      void merge_done(std::size_t p)
      {
        std::unique_lock<spinlock> _(m_lock_1);
        m_merge[p] = 1;
      }

    private:
      bool m_config {false};
      bool m_repart {false};
      std::vector<std::uint8_t> m_superk;
      std::vector<std::uint8_t> m_count;
      std::vector<std::uint8_t> m_merge;
      std::size_t m_nb;
      std::size_t m_part;
      std::string m_path;

      spinlock m_lock_1;
      spinlock m_lock_2;
  };
}

#endif /* end of include guard: STATE_HPP_OSPLYBCL */
