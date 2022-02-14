#include <kmtricks/howde/bit_vector.hpp>

namespace km
{


/* bit_vector implem */

bit_vector::bit_vector() {}

bit_vector::bit_vector(size_t nbits)
  : m_size_bits(nbits)
{
  m_bits = std::make_unique<sdsl::bit_vector(m_size_bits, 0);
  m_is_loaded = true;
}

std::string bit_vector::type() const
{
  return "bit_vector";
}

size_t bit_vector::size() const
{
  if (m_is_loaded)
    return m_bits->size();
  return 0;
}

void bit_vector::set_one(uint64_t index)
{
  (*m_bits)[index] = 1;
}

void bit_vector::set_zero(uint64_t index)
{
  (*m_bits)[index] = 0;
}

void bit_vector::set_bit(uint64_t index, bool value)
{
  (*m_bits)[index] = value;
}

bool bit_vector::check_bit(uint64_t index) const
{
  return (*m_bits)[index];
}

uint64_t bit_vector::rank1(uint64_t index) const
{
  return m_rank1->rank(index);
}

uint64_t bit_vector::select0(uint64_t rank) const
{
  return m_select0->select(rank + 1);
}

void bit_vector::serialize(std::ostream& out_stream)
{
  if (out_stream.good())
    m_bits->serialize(out_stream, nullptr, "");
  else
    std::cout << "ok" << std::endl;
}

void bit_vector::deserialize(std::ostream& in_stream)
{
  if (!m_is_loaded)
  {
    m_bits = std::make_unique<sdsl::bit_vector>();
    sdsl::load(*m_bits, in_stream);
    m_size_bits = m_bits->size();
    m_is_loaded = true;
  }
}

/* rrr_vector implem */

rrr_vector::rrr_vector() {}

rrr_vector::rrr_vector(size_t nbits)
  : m_size_bits(nbits)
{
  m_bits = std::make_unique<sdsl::bit_vector>(m_size_bits, 0);
  m_is_loaded = true;
}

std::string rrr_vector::type() const
{
  return "rrr_vector";
}

size_t rrr_vector::size() const
{
  if (m_rrr)
    return m_rrr->size();
  if (m_bits)
    return m_bits->size();
  return 0;
}

void rrr_vector::set_one(uint64_t index)
{
  base::set_one(index);
}

void rrr_vector::set_zero(uint64_t index)
{
  base::set_zero(index);
}

void rrr_vector::set_bit(uint64_t index, bool value)
{
  base::set_bit(index, value);
}

bool rrr_vector::check_bit(uint64_t index) const
{
  if (m_rrr)
    return (*m_rrr)[index];
  else
    return base::check_bit(index);
}

uint64_t rrr_vector::rank1(uint64_t index) const
{
  return m_rrr_rank1->rank(index);
}

uint64_t rrr_vector::select0(uint64_t rank) const
{
  return m_rrr_select0->select(rank + 1);
}

void rrr_vector::serialize(std::ostream& out_stream)
{
  if (!m_finished)
  {
    if (m_bits)
      base::serialize(out_stream);
  }
  else
  {
    if (!m_rrr)
      compress();
    m_rrr->serialize(out_stream, nullptr, "");
  }
}

void rrr_vector::deserialize(std::ostream& in_stream)
{

}

}
