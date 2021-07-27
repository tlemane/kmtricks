#pragma once
#include <string>
#include <vector>
#include <kmtricks/kmer.hpp>

namespace km {

template<size_t MAX_K>
class SuperKmer
{
 inline static size_t s_kmer_size;
public:
  SuperKmer() {}

  SuperKmer(uint8_t* buffer, size_t size, size_t kmer_size)
  {
    set(buffer, size, kmer_size);
  }

  SuperKmer(const std::string& superk, size_t kmer_size)
  {
    set(superk, kmer_size);
  }

  void set_k(size_t k)
  {
    s_kmer_size = k;
  }

  void set_size(size_t size)
  {
    m_size = size;
  }

  void set(const std::string& superk, size_t kmer_size)
  {
    set_k(kmer_size);
    set_size(superk.size());
    if ((m_data.size() / 4) < m_size)
      m_data.resize((m_size / 4) + 1, 0);
    for (size_t i=0; i<m_size; i++)
    {
      int pos = i/4;
      m_data[pos] <<= 2;
      m_data[pos] |= NToB[superk[i]];
    }
    int shift = m_size % 4;
    m_data[(m_size - 1)/4] <<= shift ? 8 - (m_size % 4) * 2 : 0;
  }

  void set(uint8_t* buffer, size_t size, size_t k)
  {
    set_k(k); set_size(size);
    m_data.resize(0, (size / 4) + 1);
    uint8_t* ptr = buffer;
    uint8_t nbK = size - s_kmer_size + 1;
    uint8_t newbyte;
    int rem_size = s_kmer_size;

    uint64_t Tnewbyte;
    int nbr = 0;
    int index = s_kmer_size % 4 ? (s_kmer_size / 4) : (s_kmer_size / 4) - 1;
    while (rem_size >= 4)
    {
      newbyte = *ptr; ptr++;
      m_data[index] = newbyte; index--;
      Tnewbyte = newbyte;
      rem_size -= 4; nbr++;
    }

    int uid = 4;

    if (rem_size > 0)
    {
      newbyte = *ptr; ptr++;
      m_data[index] = newbyte;
      Tnewbyte = newbyte;
      uid = rem_size;
    }

    uint8_t rem = nbK;
    size_t c_offset = s_kmer_size / 4;
    size_t shift = s_kmer_size % 4 ? (6-(2*((s_kmer_size%4)-1))) : 6;

    m_data[c_offset] >>= shift;

    int nbnt = s_kmer_size % 4;

    for (int i=0; i<nbK; i++, rem--)
    {
      if (rem < 2) break;
      if (uid >= 4)
      {
        newbyte = *ptr; ptr++;
        Tnewbyte = newbyte;
        uid = 0;
      }
      uint8_t newnt = (Tnewbyte >> (2*uid)) & 3; uid++;
      m_data[c_offset] <<= 2;
      m_data[c_offset] |= newnt;
      nbnt++;
      if (nbnt == 4)
      {
        c_offset++;
        nbnt = 0;
      }
    }

    while (nbnt < 4)
    {
      m_data[c_offset] <<= 2;
      nbnt++;
    }
  }

  std::string to_string()
  {
    std::string str;
    size_t nb_bytes = (m_size/4) + 1;
    for (size_t i=0; i<nb_bytes; i++)
    {
      uint8_t tmp = m_data[i];
      for (int j=3; j>=0; j--)
      {
        str += bToN[(tmp >> (2 * j)) & 3];
      }
    }
    return str.substr(0, m_size);
  }

private:
  std::vector<uint8_t> m_data;
  size_t m_size;
};

};