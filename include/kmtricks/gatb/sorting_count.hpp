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
#include <gatb/gatb_core.hpp>
#include <gatb/bank/api/IBank.hpp>
#include <gatb/kmer/api/ICountProcessor.hpp>
#include <gatb/tools/designpattern/api/ICommand.hpp>
#include <gatb/tools/collections/impl/Hash16.hpp>
#include <gatb/tools/misc/api/Abundance.hpp>
#include <robin_hood.h>

#include <kmtricks/gatb/count_processor.hpp>
#include <kmtricks/superk.hpp>

#include <sabuhash.h>

#include <spdlog/spdlog.h>

#include <xxhash.h>

#define IX(x, rad) ((rad)+(256)*(x))

namespace km {

class CounterBuilder
{
public:
  CounterBuilder(size_t nbBanks = 1) : _abundancePerBank(nbBanks) {}

  size_t size() const { return _abundancePerBank.size(); }

  void init(size_t idxBank = 0)
  {
    for (size_t k = 0; k < _abundancePerBank.size(); k++)
    {
      _abundancePerBank[k] = 0;
    }
    _abundancePerBank[idxBank] = 1;
  }

  void increase(size_t idxBank = 0) { _abundancePerBank[idxBank]++; }

  void set(CountNumber val, size_t idxBank = 0) { _abundancePerBank[idxBank] = val; }

  CountNumber operator[](size_t idxBank) const { return _abundancePerBank[idxBank]; }

  const CountVector &get() const { return _abundancePerBank; }

private:
  CountVector _abundancePerBank;
};


template <typename CountProcessor, typename Storage, size_t span>
class IPartitionCounter
{
  typedef typename ::Kmer<span>::Type Type;
  typedef typename ::Kmer<span>::Count Count;

public:
  IPartitionCounter(CountProcessor *processor,
                    size_t kmer_size,
                    PartiInfo<5>* pinfo,
                    MemAllocator& pool,
                    Storage *superk_storage,
                    uint32_t part)
      : m_processor(processor), m_kmer_size(kmer_size), m_pinfo(pinfo), m_pool(pool),
        m_superk_storage(superk_storage), m_part(part)
  {
  }

protected:
  void insert(const Type &kmer, const CounterBuilder &count)
  {
    m_processor->process(m_part, kmer, count.get());
  }

  void insert(const Type&kmer, uint32_t count)
  {
    m_processor->process(m_part, kmer, count);
  }

  void insert_hash(uint64_t hash, const CounterBuilder &count)
  {
    m_processor->process(m_part, hash, count.get());
  }

  void insert_hash(uint64_t hash, uint32_t count)
  {
    m_processor->process(m_part, hash, count);
  }

  void set_processor(CountProcessor *processor)
  {
    m_processor = processor;
  }

protected:
  CountProcessor *m_processor;
  size_t m_kmer_size;
  PartiInfo<5>* m_pinfo;
  MemAllocator& m_pool;
  Storage *m_superk_storage;
  uint32_t m_part;
};

template <typename Storage, size_t span>
class ReadSuperk
{
  typedef typename ::Kmer<span>::Type Type;

public:
  ReadSuperk(Storage *superk_storage, int file_id, size_t kmer_size, uint64_t *r_idx,
             Type **radix_kmers, uint64_t *radix_sizes)
      : m_superk_storage(superk_storage), m_file_id(file_id), m_kmer_size(kmer_size),
        m_r_idx(r_idx), m_radix_kmers(radix_kmers), m_radix_sizes(radix_sizes)
  {
    m_kx = 4;
    Type un;
    un.setVal(1);
    m_kmer_mask = (un << (m_kmer_size * 2)) - 1;
    m_mask_radix.setVal((int64_t)255);
    m_mask_radix = m_mask_radix << ((m_kmer_size - 4) * 2);
    m_shift = 2 * (m_kmer_size - 1);
    m_shift_val = un.getSize() - 8;
    m_shift_radix = ((m_kmer_size - 4) * 2);
  }

  void execute()
  {
    uint nbbreak = 0;
    unsigned int nb_bytes_read;
    while (m_superk_storage->readBlock(&m_buffer, &m_buffer_size, &nb_bytes_read, m_file_id))
    {
      //decode block and iterate through its superkmers
      unsigned char *ptr = m_buffer;
      u_int8_t nbK; //number of kmers in the superkmer
      int nbsuperkmer_read = 0;
      u_int8_t newbyte = 0;

      while (ptr < (m_buffer + nb_bytes_read)) //decode whole block
      {
        //decode a superkmer
        nbK = *ptr;
        ptr++;
        //int nb_bytes_superk = (_kmerSize + nbK -1 +3) /4  ;

        int rem_size = m_kmer_size;

        Type Tnewbyte;
        int nbr = 0;
        m_seedk.setVal(0);
        while (rem_size >= 4)
        {
          newbyte = *ptr;
          ptr++;
          Tnewbyte.setVal(newbyte);

          m_seedk = m_seedk | (Tnewbyte << (8 * nbr));
          rem_size -= 4;
          nbr++;
        }

        int uid = 4; //uid = nb nt used in current newbyte

        //reste du seed kmer
        if (rem_size > 0)
        {
          newbyte = *ptr;
          ptr++;
          Tnewbyte.setVal(newbyte);

          m_seedk = (m_seedk | (Tnewbyte << (8 * nbr)));
          uid = rem_size;
        }
        m_seedk = m_seedk & m_kmer_mask;

        ///////////////////////// seedk should be ready here , now parse kx-mers ////////////////////////////////

        u_int8_t rem = nbK;
        Type temp = m_seedk;
        Type rev_temp = revcomp(temp, m_kmer_size);
        Type newnt;
        Type mink, prev_mink;
        prev_mink.setVal(0);
        uint64_t idx;

        bool prev_which = (temp < rev_temp);

        int kx_size = -1; //next loop start at ii=0, first kmer will put it at 0
        Type radix_kxmer_forward = (temp & m_mask_radix) >> ((m_kmer_size - 4) * 2);
        Type first_revk, kinsert, radix_kxmer;
        first_revk.setVal(0);

        if (!prev_which)
          first_revk = rev_temp;

        u_int8_t rid;

        for (int ii = 0; ii < nbK; ii++, rem--)
        {
          bool which = (temp < rev_temp);
          mink = which ? temp : rev_temp;

          nbbreak++;
          if (which != prev_which || kx_size >= m_kx) // kxmer_size = 1
          {
            //output kxmer size kx_size,radix_kxmer
            //kx mer is composed of superKp[ii-1] superKp[ii-2] .. superKp[ii-n] with nb elems  n  == kxmer_size +1  (un seul kmer ==k+0)

            if (prev_which)
            {
              radix_kxmer = radix_kxmer_forward;
              kinsert = prev_mink;
            }
            else // si revcomp, le radix du kxmer est le debut du dernier kmer
            {
              //previous mink
              radix_kxmer = (prev_mink & m_mask_radix) >> m_shift_radix;
              kinsert = first_revk;
            }

            //record kxmer
            rid = radix_kxmer.getVal();
            //idx = _r_idx[IX(kx_size,rid)]++;
            idx = __sync_fetch_and_add(m_r_idx + IX(kx_size, rid), 1); // si le sync fetch est couteux, faire un mini buffer par thread

            m_radix_kmers[IX(kx_size, rid)][idx] = kinsert << ((4 - kx_size) * 2); //[kx_size][rid]

            radix_kxmer_forward = (mink & m_mask_radix) >> m_shift_radix;
            kx_size = 0;

            if (!which)
              first_revk = rev_temp;
          }
          else
          {
            kx_size++;
          }

          prev_which = which;
          prev_mink = mink;

          if (rem < 2)
            break; //no more kmers in this superkmer, the last one has just been eaten

          //////////////////////////////now decode next kmer of this superkmer //////////////////////////////////////////////

          if (uid >= 4) //read next byte
          {
            newbyte = *ptr;
            ptr++;
            Tnewbyte.setVal(newbyte);
            uid = 0;
          }

          newnt = (Tnewbyte >> (2 * uid)) & 3;
          uid++;

          temp = ((temp << 2) | newnt) & m_kmer_mask;

          newnt.setVal(comp_NT[newnt.getVal()]);
          rev_temp = ((rev_temp >> 2) | (newnt << m_shift)) & m_kmer_mask;

          ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }

        //record last kxmer prev_mink et monk ?
        if (prev_which)
        {
          radix_kxmer = radix_kxmer_forward;
          kinsert = prev_mink;
        }
        else // si revcomp, le radix du kxmer est le debut du dernier kmer
        {
          //previous mink
          radix_kxmer = (prev_mink & m_mask_radix) >> m_shift_radix;
          kinsert = first_revk;
        }

        //record kxmer
        rid = radix_kxmer.getVal();
        //idx = _r_idx[IX(kx_size,rid)]++;
        idx = __sync_fetch_and_add(m_r_idx + IX(kx_size, rid), 1); // si le sync fetch est couteux, faire un mini buffer par thread

        m_radix_kmers[IX(kx_size, rid)][idx] = kinsert << ((4 - kx_size) * 2); // [kx_size][rid]
        //cout << "went okay " << idx << endl;

        //////////////////////////////////////////////////////////
        //ptr+=nb_bytes_superk;

        //now go to next superk of this block, ptr should point to beginning of next superk
        nbsuperkmer_read++;
        /////////
      }
    }

    if (m_buffer != 0)
      free(m_buffer);
  }

private:
  Storage* m_superk_storage;
  int m_file_id;
  unsigned char *m_buffer {0};
  unsigned int m_buffer_size {0};
  size_t m_kmer_size;
  int m_kx;

  Type **m_radix_kmers;
  uint64_t *m_radix_sizes;
  uint64_t *m_r_idx;

  Type m_superk, m_seedk;
  Type m_radix, m_mask_radix;
  Type m_kmer_mask;
  size_t m_shift;
  size_t m_shift_val;
  size_t m_shift_radix;
};

template <size_t span>
struct IHasher
{
  typedef typename ::Kmer<span>::Type Type;
  static constexpr size_t slot {(span + 31) / 32};
  virtual uint64_t operator()(Type& kmer) = 0;
  virtual ~IHasher() = default;
};

template<size_t span>
using hasher_t = std::unique_ptr<IHasher<span>>;

template <size_t span>
struct KmSabuhash : public IHasher<span>
{
public:
  KmSabuhash(size_t kmer_size, uint64_t win, uint64_t p)
    : m_kmer_size(kmer_size), m_win(win), m_p(p), m_hasher(SabuHash(m_kmer_size))
  {}

  typedef typename ::Kmer<span>::Type Type;
  uint64_t operator()(Type &kmer)
  {
    return (m_hasher.hash(kmer.toString(m_kmer_size)) % m_win) + (m_win * m_p);
  }
private:
  size_t m_kmer_size;
  SabuHash m_hasher;
  uint64_t m_win;
  uint64_t m_p;
};

template <size_t span>
struct KmXXHash : public IHasher<span>
{
public:
  KmXXHash(size_t kmer_size, uint64_t win, uint64_t p)
    : m_kmer_size(kmer_size), m_win(win), m_p(p), m_len(((kmer_size+31)/32)*8)
  {}
  typedef typename ::Kmer<span>::Type Type;
  uint64_t operator()(Type &kmer)
  {
    return (XXH64(kmer.get_data(), m_len, 0) % m_win) + (m_win * m_p);
  }
private:
  size_t m_kmer_size;
  uint64_t m_win;
  uint64_t m_p;
  uint64_t m_len;
};

template <typename Storage, size_t span>
class ReadSuperkHash
{
  typedef typename ::Kmer<span>::Type Type;

public:
  ReadSuperkHash(Storage *superk_storage,
                 int file_id,
                 int kmer_size,
                 uint64_t *r_idx,
                 uint64_t *array,
                 uint64_t window)
      : superk_storage(superk_storage), file_id(file_id), buffer(0), buffer_size(0),
        kmer_size(kmer_size), r_idx(r_idx), array(array), win_size(window)
  {
    Type un;
    un.setVal(1);
    kmer_mask = (un << (kmer_size * 2)) - 1;
    shift = 2 * (kmer_size - 1);
    hasher = std::make_unique<KmXXHash<span>>(kmer_size, win_size, file_id);
  }

  void execute()
  {
    uint nbbreak = 0;
    uint32_t nb_bytes_read;
    while (superk_storage->readBlock(&buffer, &buffer_size, &nb_bytes_read, file_id))
    {
      unsigned char *ptr = buffer;
      uint8_t nbK;
      int nbsuperkmer_read = 0;
      uint8_t newbyte = 0;
      while (ptr < (buffer + nb_bytes_read))
      {
        nbK = *ptr;
        ptr++;
        int rem_size = kmer_size;
        Type Tnewbyte;
        int nbr = 0;
        seedk.setVal(0);
        while (rem_size >= 4)
        {
          newbyte = *ptr;
          ptr++;
          Tnewbyte.setVal(newbyte);
          seedk = seedk | (Tnewbyte << (8 * nbr));
          rem_size -= 4;
          nbr++;
        }
        int uid = 4;
        if (rem_size > 0)
        {
          newbyte = *ptr;
          ptr++;
          Tnewbyte.setVal(newbyte);
          seedk = (seedk | (Tnewbyte << (8 * nbr)));
          uid = rem_size;
        }
        seedk = seedk & kmer_mask;

        uint8_t rem = nbK;
        Type temp = seedk;
        Type rev_temp = revcomp(temp, kmer_size);
        Type mink, newnt;
        uint64_t idx;

        bool which = (temp < rev_temp);
        mink = which ? temp : rev_temp;

        idx = *r_idx;
        array[idx] = (*hasher.get())(mink);
        ++*r_idx;

        for (int i = 0; i < nbK; i++, rem--)
        {
          nbbreak++;
          if (rem < 2)
            break;

          if (uid >= 4)
          {
            newbyte = *ptr;
            ptr++;
            Tnewbyte.setVal(newbyte);
            uid = 0;
          }

          newnt = (Tnewbyte >> (2 * uid)) & 3;
          uid++;
          temp = ((temp << 2) | newnt) & kmer_mask;
          newnt.setVal(comp_NT[newnt.getVal()]);
          rev_temp = ((rev_temp >> 2) | (newnt << shift)) & kmer_mask;

          which = (temp < rev_temp);
          mink = which ? temp : rev_temp;
          idx = *r_idx;

          array[idx] = (*hasher.get())(mink);
          ++*r_idx;
        }
        nbsuperkmer_read++;
      }
    }
    if (buffer != 0)
      free(buffer);
  }

private:
  Storage *superk_storage;
  int file_id;
  unsigned char *buffer;
  unsigned int buffer_size;
  int kmer_size;

  uint64_t *r_idx;
  Type seedk;
  Type kmer_mask;
  size_t shift;
  uint64_t *array;
  uint64_t win_size;
  hasher_t<span> hasher;
};

template <size_t span>
class KmerSort
{
public:
  typedef typename ::Kmer<span>::Type Type;
  KmerSort(Type **kmer_vector, int begin, int end, uint64_t* radix_size)
    : begin(begin), end(end), m_kmer_vector(kmer_vector), m_radix_size(radix_size)
  {
  }

  void execute()
  {
    for (int ii = begin; ii <= end; ii++)
    {
      if (m_radix_size[ii] > 0)
      {
        Type *kmers = m_kmer_vector[ii];
        std::sort(&kmers[0], &kmers[m_radix_size[ii]]);
      }
    }
  }

private:
  int begin;
  int end;
  Type **m_kmer_vector;
  uint64_t* m_radix_size;
};

class HashSort
{
public:
  HashSort(uint64_t* hash_vector, size_t array_size)
      : m_hash_vector(hash_vector), m_size(array_size)
  {
  }

  void execute()
  {
    std::sort(&m_hash_vector[0], &m_hash_vector[m_size]);
  }

private:
  uint64_t* m_hash_vector;
  size_t m_size;
};

template <size_t span>
class KXmerPointer
{
  typedef typename ::Kmer<span>::Type Type;

public:
  KXmerPointer(Type **kxmers,
               int prefix_size,
               int x_size,
               int min_radix,
               int max_radix,
               int kmer_size,
               uint64_t *radix_sizes)
      : m_kxmers(kxmers), m_radix_sizes(radix_sizes), m_prefix_size(prefix_size), m_x_size(x_size),
        m_low_radix(min_radix), m_high_radix(max_radix), m_kmer_size(kmer_size),
        m_idx_radix(min_radix)
  {
    Type un;
    un.setVal(1);
    m_kmer_mask = (un << (m_kmer_size * 2)) - un;
    m_shift_size = ((4 - m_prefix_size) * 2);
    m_radix_mask.setVal(m_idx_radix);
    m_radix_mask = m_radix_mask << ((m_kmer_size - 4) * 2);
    m_radix_mask = m_radix_mask << (2 * m_prefix_size);
  }

  bool next()
  {
    m_cur_idx++;
    while (m_idx_radix <= m_high_radix &&
           static_cast<uint64_t>(m_cur_idx) >= m_radix_sizes[m_idx_radix])
    {
      m_idx_radix++;
      m_cur_idx = 0;
      m_radix_mask.setVal(m_idx_radix);
      m_radix_mask = m_radix_mask << ((m_kmer_size - 4) * 2);
      m_radix_mask = m_radix_mask << (2 * m_prefix_size);
    }
    return (m_idx_radix <= m_high_radix);
  }

  Type value() const
  {
    return (((m_kxmers[m_idx_radix][m_cur_idx]) >> m_shift_size) | m_radix_mask) & m_kmer_mask;
  }

  uint32_t getBankId() const
  {
    return 0;
  }

private:
  Type **m_kxmers;
  uint64_t *m_radix_sizes;
  int64_t m_cur_idx{-1};
  Type m_cur_radix;
  Type m_kmer_mask;
  Type m_radix_mask;
  int m_idx_radix;
  int m_low_radix;
  int m_high_radix;
  int m_shift_size;
  int m_prefix_size;
  int m_kmer_size;
  int m_x_size;
};

template <typename Storage, size_t span>
class KmerPartCounter : public IPartitionCounter<ICountProcessor<span>, Storage, span>
{
public:
  typedef typename ::Kmer<span>::Type Type;
  typedef typename ::Kmer<span>::Count Count;
  typedef ICountProcessor<span> CountProcessor;
  static const size_t KX = 4;

private:
  typedef std::pair<int, Type> kxp;
  struct kxpcomp
  {
    bool operator()(kxp l, kxp r)
    {
      return r.second < l.second;
    }
  };

public:
  KmerPartCounter(CountProcessor *processor,
                  PartiInfo<5>* pinfo,
                  int parti,
                  size_t kmer_size,
                  MemAllocator &pool,
                  Storage *superk_storage)
      : IPartitionCounter<CountProcessor, Storage, span>(processor,
                                                kmer_size,
                                                pinfo, pool,
                                                superk_storage,
                                                parti),
        radix_kmers(0), radix_sizes(0), r_idx(0)
  {
  }

  void execute()
  {
    radix_kmers = (Type **)MALLOC(256 * (KX + 1) * sizeof(Type *));
    radix_sizes = (uint64_t *)MALLOC(256 * (KX + 1) * sizeof(uint64_t));
    r_idx = (uint64_t *)CALLOC(256 * (KX + 1), sizeof(uint64_t));

    executeRead();
    executeSort();
    executeDump();

    FREE(radix_kmers);
    FREE(radix_sizes);
    FREE(r_idx);
  }

private:
  void executeRead()
  {
    if constexpr(std::is_same_v<Storage, SuperKmerBinFiles>)
      this->m_superk_storage->openFile("r", this->m_part);
    else
      this->m_superk_storage->openFile(this->m_part);

    uint64_t sum_nbxmer = 0;

    {
      LocalSynchronizer synchro(this->m_pool.getSynchro());
      this->m_pool.align(16);

      for (size_t xx = 0; xx < (KX + 1); xx++)
      {
        for (int ii = 0; ii < 256; ii++)
        {
          size_t nb_kmers = this->m_pinfo->getNbKmer(this->m_part, ii, xx);
          radix_kmers[IX(xx, ii)] = (Type *)this->m_pool.pool_malloc(nb_kmers * sizeof(Type),
                                                                     "kmers alloc");
          radix_sizes[IX(xx, ii)] = nb_kmers;
          sum_nbxmer += nb_kmers;
        }
      }

      ReadSuperk<Storage, span> read_cmd(this->m_superk_storage, this->m_part, this->m_kmer_size,
                                r_idx, radix_kmers, radix_sizes);
      read_cmd.execute();
    }
    this->m_superk_storage->closeFile(this->m_part);
  }

  void executeSort()
  {
    for (size_t xx = 0; xx < (KX + 1); xx++)
    {
      KmerSort<span> sort_cmd(radix_kmers + IX(xx, 0), 0, 255, radix_sizes + IX(xx, 0));
      sort_cmd.execute();
    }
  }

  void executeDump()
  {
    int nbkxpointers = 453;
    std::vector<KXmerPointer<span> *> vec_pointer(nbkxpointers);
    int best_p;
    std::priority_queue<kxp, std::vector<kxp>, kxpcomp> pq;
    size_t nbBanks = nb_items_per_bank_per_part.size();

    if (nbBanks == 0)
      nbBanks = 1;

    CounterBuilder solidCounter(nbBanks);
    Type previous_kmer;
    int pidx = 0;

    vec_pointer[pidx++] = new KXmerPointer<span>(radix_kmers + IX(0, 0), 0, 0, 0, 255,
                                                 this->m_kmer_size, radix_sizes + IX(0, 0));

    vec_pointer[pidx++] = new KXmerPointer<span>(radix_kmers + IX(1, 0), 0, 1, 0, 255,
                                                 this->m_kmer_size, radix_sizes + IX(1, 0));
    int lowr = 0;
    int maxr = 63;

    for (uint32_t ii = 0; ii < 4; ii++)
    {
      vec_pointer[pidx++] = new KXmerPointer<span>(radix_kmers + IX(1, 0), 1, 1, lowr, maxr,
                                                   this->m_kmer_size, radix_sizes + IX(1, 0));
      lowr += 64;
      maxr += 64;
    }

    vec_pointer[pidx++] = new KXmerPointer<span>(radix_kmers + IX(2, 0), 0, 2, 0, 255,
                                                 this->m_kmer_size, radix_sizes + IX(2, 0));
    lowr = 0;
    maxr = 63;
    for (uint32_t ii = 0; ii < 4; ii++)
    {
      vec_pointer[pidx++] = new KXmerPointer<span>(radix_kmers + IX(2, 0), 1, 2, lowr, maxr,
                                                   this->m_kmer_size, radix_sizes + IX(2, 0));
      lowr += 64;
      maxr += 64;
    }

    lowr = 0;
    maxr = 15;
    for (uint32_t ii = 0; ii < 16; ii++)
    {
      vec_pointer[pidx++] = new KXmerPointer<span>(radix_kmers + IX(2, 0), 2, 2, lowr, maxr,
                                                   this->m_kmer_size, radix_sizes + IX(2, 0));
      lowr += 16;
      maxr += 16;
    }

    vec_pointer[pidx++] = new KXmerPointer<span>(radix_kmers + IX(3, 0), 0, 3, 0, 255,
                                                 this->m_kmer_size, radix_sizes + IX(3, 0));

    lowr = 0;
    maxr = 63;
    for (uint32_t ii = 0; ii < 4; ii++)
    {
      vec_pointer[pidx++] = new KXmerPointer<span>(radix_kmers + IX(3, 0), 1, 3, lowr, maxr,
                                                   this->m_kmer_size, radix_sizes + IX(3, 0));
      lowr += 64;
      maxr += 64;
    }

    lowr = 0;
    maxr = 15;
    for (uint32_t ii = 0; ii < 16; ii++)
    {
      vec_pointer[pidx++] = new KXmerPointer<span>(radix_kmers + IX(3, 0), 2, 3, lowr, maxr,
                                                   this->m_kmer_size, radix_sizes + IX(3, 0));
      lowr += 16;
      maxr += 16;
    }

    lowr = 0;
    maxr = 3;
    for (uint32_t ii = 0; ii < 64; ii++)
    {
      vec_pointer[pidx++] = new KXmerPointer<span>(radix_kmers + IX(3, 0), 3, 3, lowr, maxr,
                                                   this->m_kmer_size, radix_sizes + IX(3, 0));
      lowr += 4;
      maxr += 4;
    }

    vec_pointer[pidx++] = new KXmerPointer<span>(radix_kmers + IX(4, 0), 0, 4, 0, 255,
                                                 this->m_kmer_size, radix_sizes + IX(4, 0));

    lowr = 0;
    maxr = 63;
    for (uint32_t ii = 0; ii < 4; ii++)
    {
      vec_pointer[pidx++] = new KXmerPointer<span>(radix_kmers + IX(4, 0), 1, 4, lowr, maxr,
                                                   this->m_kmer_size, radix_sizes + IX(4, 0));
      lowr += 64;
      maxr += 64;
    }

    lowr = 0;
    maxr = 15;
    for (uint32_t ii = 0; ii < 16; ii++)
    {
      vec_pointer[pidx++] = new KXmerPointer<span>(radix_kmers + IX(4, 0), 2, 4, lowr, maxr,
                                                   this->m_kmer_size, radix_sizes + IX(4, 0));
      lowr += 16;
      maxr += 16;
    }

    lowr = 0;
    maxr = 3;
    for (uint32_t ii = 0; ii < 64; ii++)
    {
      vec_pointer[pidx++] = new KXmerPointer<span>(radix_kmers + IX(4, 0), 3, 4, lowr, maxr,
                                                   this->m_kmer_size, radix_sizes + IX(4, 0));
      lowr += 4;
      maxr += 4;
    }

    lowr = 0;
    maxr = 0;
    for (uint32_t ii = 0; ii < 256; ii++)
    {
      vec_pointer[pidx++] = new KXmerPointer<span>(radix_kmers + IX(4, 0), 4, 4, lowr, maxr,
                                                   this->m_kmer_size, radix_sizes + IX(4, 0));
      lowr += 1;
      maxr += 1;
    }

    for (int ii = 0; ii < nbkxpointers; ii++)
    {
      if (vec_pointer[ii]->next())
      {
        pq.push(kxp(ii, vec_pointer[ii]->value()));
      }
    }

    uint32_t count = 0;
    if (pq.size() != 0)
    {
      best_p = pq.top().first;
      pq.pop();
      previous_kmer = vec_pointer[best_p]->value();
      //solidCounter.init(vec_pointer[best_p]->getBankId());
      count = 1;
      while (1)
      {
        if (!vec_pointer[best_p]->next())
        {
          if (pq.size() == 0)
            break;
          best_p = pq.top().first;
          pq.pop();
        }

        if (vec_pointer[best_p]->value() != previous_kmer)
        {
          pq.push(kxp(best_p, vec_pointer[best_p]->value()));
          best_p = pq.top().first;
          pq.pop();

          if (vec_pointer[best_p]->value() != previous_kmer)
          {
            //this->insert(previous_kmer, solidCounter);
            this->insert(previous_kmer, count);
            //solidCounter.init(vec_pointer[best_p]->getBankId());
            count = 1;
            previous_kmer = vec_pointer[best_p]->value();
          }
          else
          {
            //solidCounter.increase(vec_pointer[best_p]->getBankId());
            count++;
          }
        }
        else
        {
          //solidCounter.increase(vec_pointer[best_p]->getBankId());
          count++;
          //solidCounter.increase(vec_pointer[best_p]->getBankId());
        }
      }
      //this->insert(previous_kmer, solidCounter);
      this->insert(previous_kmer, count);
    }

    for (int ii = 0; ii < nbkxpointers; ii++)
    {
      delete vec_pointer[ii];
    }
  }

private:
  Type **radix_kmers;
  uint64_t *radix_sizes;
  uint64_t *r_idx;
  std::vector<size_t> nb_items_per_bank_per_part;
};

//template <size_t span>
//class KmerPartCounterByHash : public IPartitionCounter<ICountProcessor<span>, span>
//{
//public:
//  typedef typename ::Kmer<span>::Type Type;
//  typedef typename ::Kmer<span>::Count Count;
//  typedef ICountProcessor<span> CountProcessor;
//
//private:
//  uint64_t m_hash_memory;
//  uint64_t m_window;
//  uint64_t m_part_id;
//};


template <typename Storage, size_t span>
class HashPartCounter : public IPartitionCounter<IHashProcessor<span>, Storage, span>
{
public:
  typedef typename ::Kmer<span>::Type Type;
  typedef typename ::Kmer<span>::Count Count;
  typedef IHashProcessor<span> CountProcessor;
  static const size_t KX = 4;

public:
  HashPartCounter(CountProcessor *processor,
                  PartiInfo<5>* pinfo,
                  int parti,
                  size_t kmer_size,
                  MemAllocator &pool,
                  Storage *superk_storage,
                  uint64_t window)
      : IPartitionCounter<CountProcessor, Storage, span>(processor,
                                                kmer_size,
                                                pinfo,
                                                pool,
                                                superk_storage,
                                                parti), r_idx(0), window(window)
  {
  }

  void execute()
  {
    r_idx = (uint64_t *)CALLOC(256 * (KX + 1), sizeof(uint64_t));

    executeRead();
    executeSort();
    executeDump();
    this->m_processor->finish();
    FREE(r_idx);
  }

private:
  void executeRead()
  {
    if constexpr(std::is_same_v<Storage, SuperKmerBinFiles>)
      this->m_superk_storage->openFile("r", this->m_part);
    else
      this->m_superk_storage->openFile(this->m_part);

    LocalSynchronizer synchro(this->m_pool.getSynchro());
    this->m_pool.align(16);

    size_t nb_kmers = this->m_pinfo->getNbKmer(this->m_part);
    array = (uint64_t*) this->m_pool.pool_malloc(nb_kmers*sizeof(uint64_t), std::to_string(this->m_part).c_str());
    ReadSuperkHash<Storage, span> read_cmd(this->m_superk_storage, this->m_part,
                                  this->m_kmer_size, r_idx, array, window);
    read_cmd.execute();

    this->m_superk_storage->closeFile(this->m_part);
  }

  void executeSort()
  {
    HashSort sort_cmd(array, *r_idx);
    sort_cmd.execute();
  }

  void executeDump()
  {
    uint32_t count = 1;
    uint64_t previous_kmer;
    previous_kmer = array[0];
    for (uint64_t i=1; i<*r_idx; i++)
    {
      if (previous_kmer != array[i])
      {
        this->insert_hash(previous_kmer, count);
        previous_kmer = array[i];
        count = 1;
      }
      else
      {
        count++;
      }
    }
    this->insert_hash(previous_kmer, count);
  }

private:
  uint64_t* r_idx;
  uint64_t* array;
  uint64_t window;
  std::vector<size_t> nb_items_per_bank_per_part;
};


template<typename count_type = uint16_t> struct AbundanceH
{
  AbundanceH (const uint64_t val, const uint16_t ab) : value(val), abundance(ab) {}
  AbundanceH (const uint64_t val) : value(val), abundance(0) {}
  AbundanceH() : value(0), abundance(0) {}

  AbundanceH& operator=(const AbundanceH a)
  {
    if (&a != this)
    {
      value = a.value;
      abundance = a.abundance;
    }
    return *this;
  }

  const count_type& getAbundance() const { return abundance; }
  const uint64_t& getValue() const { return value; }

  bool operator==(const AbundanceH& other) const
  {
    return value == other.value && abundance == other.abundance;
  }

  uint64_t value;
  count_type abundance;
};

template <typename Storage, size_t span>
class HashPartCounterByHash : public IPartitionCounter<IHashProcessor<span>, Storage, span>
{
public:
  typedef typename ::Kmer<span>::Type Type;
  typedef typename ::Kmer<span>::Count Count;
  typedef IHashProcessor<span> CountProcessor;

  HashPartCounterByHash(CountProcessor* processor,
                  PartiInfo<5>* pinfo,
                  int parti,
                  size_t kmer_size,
                  MemAllocator& pool,
                  Storage* superk_storage,
                  uint64_t window)
    : IPartitionCounter<IHashProcessor<span>, Storage, span>(
      processor, kmer_size, pinfo, pool, superk_storage, parti), m_window(window)
    {

    }

  void execute()
  {
    typedef typename gatb::core::tools::collections::impl::Hash16<uint64_t>::cell cell_t;
    this->m_superk_storage->openFile(this->m_part);
    Hash16<uint64_t> hash16(m_hash_memory / (1ULL << 20));
    typedef AbundanceH<uint16_t> abundance_t;

    CounterBuilder solidCounter;

    std::vector<std::string> tmp_count_filenames;
    int ks = this->m_kmer_size;
    Type un; un.setVal(1);

    Type kmerMask = (un << (ks*2)) - un;
    size_t shift = 2*(ks-1);

    Type _seedk;

    int _fileId = this->m_part;
    unsigned char * _buffer = 0 ;
    unsigned int _buffer_size = 0;

    hasher_t<span> hasher = std::make_unique<KmXXHash<span>>(this->m_kmer_size, m_window, this->m_part);
    unsigned int nb_bytes_read;
    while (this->m_superk_storage->readBlock(&_buffer, &_buffer_size, &nb_bytes_read, _fileId))
    {
      unsigned char* ptr = _buffer;
      uint8_t nbK;
      uint8_t newbyte = 0;
      while (ptr < (_buffer + nb_bytes_read))
      {
        nbK = *ptr; ptr++;
        int rem_size = this->m_kmer_size;
        Type Tnewbyte;
        int nbr = 0;
        _seedk.setVal(0);

        while (rem_size >= 4)
        {
          newbyte = *ptr; ptr++;
          Tnewbyte.setVal(newbyte);
          _seedk = _seedk | (Tnewbyte << (8 * nbr));
          rem_size -= 4; nbr++;
        }

        int uid = 4;

        if (rem_size > 0)
        {
          newbyte = *ptr; ptr++;
          Tnewbyte.setVal(newbyte);
          _seedk = (_seedk | (Tnewbyte << (8 * nbr)));
          uid = rem_size;
        }
        _seedk = _seedk & kmerMask;

        uint8_t rem = nbK;
        Type temp = _seedk;
        Type rev_temp = revcomp(temp, this->m_kmer_size);
        Type newnt;
        Type mink;

        for (int ii=0; ii<nbK; ii++, rem--)
        {
          mink = std::min(rev_temp, temp);
          hash16.insert((*hasher.get())(mink));

          if (rem < 2) break;

          if (uid >= 4)
          {
            newbyte = *ptr; ptr++;
            Tnewbyte.setVal(newbyte);
            uid = 0;
          }

          newnt = (Tnewbyte >> (2 * uid)) & 3; uid++;
          temp = ((temp << 2) | newnt) & kmerMask;
          newnt.setVal(comp_NT[newnt.getVal()]);
          rev_temp = ((rev_temp >> 2) | (newnt << shift)) & kmerMask;
        }
      }

      if (hash16.getByteSize() > m_hash_memory)
      {
        Iterator<cell_t>* itKmerAbundancePartial = hash16.iterator(true);
        LOCAL(itKmerAbundancePartial);

        std::string fname = this->m_superk_storage->getFileName(this->m_part);
        fname += gatb::core::tools::misc::impl::Stringify::format("_subpart_%i", tmp_count_filenames.size());
        tmp_count_filenames.push_back(fname);
        BagFile<abundance_t>* bagf = new BagFile<abundance_t>(fname); LOCAL(bagf);
        Bag<abundance_t>* currentbag = new BagCache<abundance_t>(bagf, 10000); LOCAL(currentbag);

        for (itKmerAbundancePartial->first(); !itKmerAbundancePartial->isDone(); itKmerAbundancePartial->next())
        {
          cell_t& cell = itKmerAbundancePartial->item();
          currentbag->insert(abundance_t(cell.graine, cell.val));
        }

        currentbag->flush();
        hash16.clear();
      }
    }

    if (_buffer != 0)
      free(_buffer);

    Iterator<cell_t>* itKmerAbundance = hash16.iterator(true); LOCAL(itKmerAbundance);

    if (tmp_count_filenames.size() != 0)
    {
      TempCountFileMerger<span> tempCountFileMerger(10, 10);
      tmp_count_filenames = tempCountFileMerger.mergeFiles(tmp_count_filenames);
      std::vector<Iterator<abundance_t>*> tmpCountIterators;

      for (uint ii=0; ii<tmp_count_filenames.size(); ii++)
      {
        std::string fname = tmp_count_filenames[ii];
        tmpCountIterators.push_back(new IteratorFile<abundance_t>(fname));
      }
      typedef std::pair<int, uint64_t> ptcf;
      struct ptcfcomp { bool operator() (ptcf l, ptcf r) { return ((r.second) < (l.second)); }};
      std::priority_queue<ptcf, std::vector<ptcf>, ptcfcomp> pq;

      itKmerAbundance->first();
      for (uint ii=0; ii<tmpCountIterators.size(); ii++)
      {
        tmpCountIterators[ii]->first();
      }

      if (!itKmerAbundance->isDone())
      {
        pq.push(ptcf(-1, itKmerAbundance->item().graine));
      }

      for (uint ii=0; ii<tmpCountIterators.size(); ii++)
      {
        if (!tmpCountIterators[ii]->isDone())
        {
          abundance_t& ab = tmpCountIterators[ii]->item();
          pq.push(ptcf(ii, ab.value));
        }
      }

      ptcf best_elem;
      int best_p;
      int current_ab = 0;
      int previous_ab = 0;
      uint64_t current, previous;

      if (pq.size() != 0)
      {
        best_elem = pq.top(); pq.pop();
        best_p = best_elem.first;
        if (best_p == -1)
        {
          previous_ab = itKmerAbundance->item().val;
        }
        else
        {
          previous_ab = tmpCountIterators[best_p]->item().abundance;
        }

        previous = best_elem.second;

        if (best_p == -1)
        {
          itKmerAbundance->next();
          if (!itKmerAbundance->isDone())
          {
            pq.push(ptcf(-1, itKmerAbundance->item().graine));
          }
        }
        else
        {
          tmpCountIterators[best_p]->next();
          if (!tmpCountIterators[best_p]->isDone())
          {
            pq.push(ptcf(best_p, tmpCountIterators[best_p]->item().value));
          }
        }

        while (pq.size() != 0)
        {
          best_elem = pq.top(); pq.pop();
          best_p = best_elem.first;

          if (best_p == -1)
          {
            current_ab = itKmerAbundance->item().val;
          }
          else
          {
            current_ab = tmpCountIterators[best_p]->item().abundance;
          }
          current = best_elem.second;

          if (best_p == -1)
          {
            itKmerAbundance->next();
            if (!itKmerAbundance->isDone())
            {
              pq.push(ptcf(-1, itKmerAbundance->item().graine));
            }
          }
          else
          {
            tmpCountIterators[best_p]->next();
            if (!tmpCountIterators[best_p]->isDone())
            {
              pq.push(ptcf(best_p, tmpCountIterators[best_p]->item().value));
            }
          }

          if (current != previous)
          {
            solidCounter.set(previous_ab);
            this->insert_hash(previous, solidCounter.get()[0]);
            previous = current;
            previous_ab = current_ab;
          }
          else
          {
            previous_ab += current_ab;
          }
        }
        solidCounter.set(previous_ab);
        this->insert_hash(previous, solidCounter.get()[0]);
      }

      for (uint ii=0; ii<tmpCountIterators.size(); ii++)
      {
        delete tmpCountIterators[ii];
      }
      for (uint ii=0; ii<tmp_count_filenames.size(); ii++)
      {
        std::string fname = tmp_count_filenames[ii];
        gatb::core::system::impl::System::file().remove(fname);
      }
    }
    else
    {
      for (itKmerAbundance->first(); !itKmerAbundance->isDone(); itKmerAbundance->next())
      {
        cell_t& cell = itKmerAbundance->item();
        solidCounter.set(cell.val);
        this->insert_hash(cell.graine, solidCounter.get()[0]);
      }
    }
    this->m_superk_storage->closeFile(this->m_part);
  }

private:
  uint64_t m_hash_memory;
  uint64_t m_window;
};

};