#include <gtest/gtest.h>
#include <kmtricks/task.hpp>
#include <kmtricks/repartition.hpp>
#include <kmtricks/gatb/gatb_utils.hpp>

#define MK 32
#define MC 4294967295

const std::string dir = "./tests_tmp/km_dir_test";
const std::string foff = "./data/kmtricks.fof";

namespace fs = std::filesystem;

int main(int argc, char* argv[])
{
  km::Kmer<MK>::m_kmer_size = 31;
  ::testing::InitGoogleTest(&argc, argv);

  fs::create_directory("./tests_tmp");
  int r = RUN_ALL_TESTS();
  fs::remove_all("./tests_tmp/km_dir_test");
  return r;
}

TEST(config_task, config_task)
{
  km::KmDir::get().init(dir, foff, true);
  {
    IProperties* props = km::get_config_properties(31, 10, 0, 0, 1, 4);
    km::ConfigTask<MK> config_task(foff, props, 5000, 4);
    config_task.exec();
  }
  {
    Storage* config_storage = StorageFactory(STORAGE_FILE).load(km::KmDir::get().m_config_storage);
    LOCAL(config_storage);

    Configuration config = Configuration();
    config.load(config_storage->getGroup("gatb"));
    ASSERT_EQ(config._kmerSize, 31);
    ASSERT_EQ(config._nb_partitions, 4);
    km::KmDir::get().init_part(config._nb_partitions);
  }
}

TEST(repart_task, repart_task)
{
  km::KmDir::get().init(dir, "", false);
  {
    km::RepartTask<MK> repart_task(foff);
    repart_task.exec();
  }

  {
    km::Repartition repart(km::KmDir::get().m_repart_storage + "_gatb/repartition.minimRepart", "");
  }
}

// Repartition depends on system configuration, so pre-computed repartition is used below
TEST(superk_task, superk_task)
{
  km::KmDir::get().init(dir, "", false);
  km::KmDir::get().m_repart_storage = "./data/repart";
  {
    std::vector<uint32_t> parts = {0, 1, 2, 3};
    km::SuperKTask<MK> superk_task("D1", true, parts);
    superk_task.exec();
    km::SuperKTask<MK> superk_task2("D2", true, parts);
    superk_task2.exec();
  }
  {
    ASSERT_TRUE(fs::exists("./tests_tmp/km_dir_test/superkmers/D1/PartiInfoFile"));
    ASSERT_TRUE(fs::exists("./tests_tmp/km_dir_test/superkmers/D1/skp.0"));
    ASSERT_TRUE(fs::exists("./tests_tmp/km_dir_test/superkmers/D1/skp.1"));
    ASSERT_TRUE(fs::exists("./tests_tmp/km_dir_test/superkmers/D1/skp.2"));
    ASSERT_TRUE(fs::exists("./tests_tmp/km_dir_test/superkmers/D1/skp.3"));
    ASSERT_TRUE(fs::exists("./tests_tmp/km_dir_test/superkmers/D1/SuperKmerBinInfoFile"));

    ASSERT_TRUE(fs::exists("./tests_tmp/km_dir_test/superkmers/D2/PartiInfoFile"));
    ASSERT_TRUE(fs::exists("./tests_tmp/km_dir_test/superkmers/D2/skp.0"));
    ASSERT_TRUE(fs::exists("./tests_tmp/km_dir_test/superkmers/D2/skp.1"));
    ASSERT_TRUE(fs::exists("./tests_tmp/km_dir_test/superkmers/D2/skp.2"));
    ASSERT_TRUE(fs::exists("./tests_tmp/km_dir_test/superkmers/D2/skp.3"));
    ASSERT_TRUE(fs::exists("./tests_tmp/km_dir_test/superkmers/D2/SuperKmerBinInfoFile"));

    {
      std::ifstream superkinfo("./tests_tmp/km_dir_test/superkmers/D1/SuperKmerBinInfoFile");
      std::string line;
      std::getline(superkinfo, line); ASSERT_EQ(line, "skp");
      std::getline(superkinfo, line);
      std::getline(superkinfo, line); ASSERT_EQ(line, "4");
      std::getline(superkinfo, line); ASSERT_EQ(line, "37");
      std::getline(superkinfo, line); ASSERT_EQ(line, "0");
      std::getline(superkinfo, line); ASSERT_EQ(line, "46");
      std::getline(superkinfo, line); ASSERT_EQ(line, "0");
      std::getline(superkinfo, line); ASSERT_EQ(line, "12");
      std::getline(superkinfo, line); ASSERT_EQ(line, "0");
      std::getline(superkinfo, line); ASSERT_EQ(line, "43");
      std::getline(superkinfo, line); ASSERT_EQ(line, "0");
    }
    {
      std::ifstream superkinfo("./tests_tmp/km_dir_test/superkmers/D2/SuperKmerBinInfoFile");
      std::string line;
      std::getline(superkinfo, line); ASSERT_EQ(line, "skp");
      std::getline(superkinfo, line);
      std::getline(superkinfo, line); ASSERT_EQ(line, "4");
      std::getline(superkinfo, line); ASSERT_EQ(line, "20");
      std::getline(superkinfo, line); ASSERT_EQ(line, "0");
      std::getline(superkinfo, line); ASSERT_EQ(line, "21");
      std::getline(superkinfo, line); ASSERT_EQ(line, "0");
      std::getline(superkinfo, line); ASSERT_EQ(line, "58");
      std::getline(superkinfo, line); ASSERT_EQ(line, "0");
      std::getline(superkinfo, line); ASSERT_EQ(line, "39");
      std::getline(superkinfo, line); ASSERT_EQ(line, "0");
    }
  }
}

TEST(count_task, kmer_count_task)
{
  km::KmDir::get().init(dir, "", false);
  Storage* config_storage = StorageFactory(STORAGE_FILE).load(km::KmDir::get().m_config_storage);
  LOCAL(config_storage);
  Configuration config = Configuration();
  config.load(config_storage->getGroup("gatb"));

  {
    km::sk_storage_t storage = std::make_shared<km::SuperKStorageReader>(km::KmDir::get().get_superk_path("D1"));
    km::parti_info_t pinfo = std::make_shared<PartiInfo<5>>(km::KmDir::get().get_superk_path("D1"));
    for (size_t p=0; p<4; p++)
    {
      std::string path = km::KmDir::get().get_count_part_path("D1", p, false, km::KM_FILE::KMER);
      km::CountTask<MK, MC, km::SuperKStorageReader> task(
        path, config, storage, pinfo, p, 0, 31, 1, false, nullptr, false);
      task.exec();
    }
  }
  {
    km::sk_storage_t storage = std::make_shared<km::SuperKStorageReader>(km::KmDir::get().get_superk_path("D2"));
    km::parti_info_t pinfo = std::make_shared<PartiInfo<5>>(km::KmDir::get().get_superk_path("D2"));
    for (size_t p=0; p<4; p++)
    {
      std::string path = km::KmDir::get().get_count_part_path("D2", p, false, km::KM_FILE::KMER);
      km::CountTask<MK, MC, km::SuperKStorageReader> task(
        path, config, storage, pinfo, p, 0, 31, 1, false, nullptr, false);
      task.exec();
    }
  }
  {
    km::Kmer<MK> kmer; kmer.set_k(31);
    uint32_t count;
    km::KmerReader<8192> kr("./tests_tmp/km_dir_test/counts/partition_0/D1.kmer");
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "AATATACTATATAATATATATAGCGAGGGGG"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ACATAATATACTATATAATATATATAGCGAG"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ACAGAGACATAATATACTATATAATATATAT"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ACAGCAGACAGAGACATAATATACTATATAA"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ACGACAGCAGACAGAGACATAATATACTATA"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ACGACGCCAGCAGAGAGACGCACACGAGACA"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ACGCCAGCAGAGAGACGCACACGAGACAGCG"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ATAATATACTATATAATATATATAGCGAGGG"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ATATATTATATAGTATATTATGTCTCTGTCT"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ATATAGCGAGGGGGGGAGAGCCAGCAGCACC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ATATAGTATATTATGTCTCTGTCTGCTGTCG"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ATATTATATAGTATATTATGTCTCTGTCTGC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ATAGCGAGGGGGGGAGAGCCAGCAGCACCCC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ATAGTATATTATGTCTCTGTCTGCTGTCGTC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ATTATATAGTATATTATGTCTCTGTCTGCTG"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "AGACATAATATACTATATAATATATATAGCG"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "AGAGACATAATATACTATATAATATATATAG"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "AGCAGACAGAGACATAATATACTATATAATA"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "AGCAGAGAGACGCACACGAGACAGCGACGAG"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "CATAATATACTATATAATATATATAGCGAGG"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "CAGACAGAGACATAATATACTATATAATATA"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "CAGAGACATAATATACTATATAATATATATA"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "CAGAGAGACGCACACGAGACAGCGACGAGCG"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "CAGCAGAGAGACGCACACGAGACAGCGACGA"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "CCAGCAGAGAGACGCACACGAGACAGCGACG"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "CCCCTCGCTATATATATTATATAGTATATTA"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "CTGTCTCGTGTGCGTCTCTCTGCTGGCGTCG"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "CGCCAGCAGAGAGACGCACACGAGACAGCGA"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "TATATATTATATAGTATATTATGTCTCTGTC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "TATATAGCGAGGGGGGGAGAGCCAGCAGCAC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "TATATAGTATATTATGTCTCTGTCTGCTGTC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "TATAGCGAGGGGGGGAGAGCCAGCAGCACCC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "TCGCTATATATATTATATAGTATATTATGTC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "GACGCCAGCAGAGAGACGCACACGAGACAGC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "GAGACATAATATACTATATAATATATATAGC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "GCAGAGAGACGCACACGAGACAGCGACGAGC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "GCCAGCAGAGAGACGCACACGAGACAGCGAC"); EXPECT_EQ(1, count);
    bool f = kr.read<MK, MC>(kmer, count); EXPECT_FALSE(f);
  }
  {
    km::Kmer<MK> kmer; kmer.set_k(31);
    uint32_t count;
    km::KmerReader<8192> kr("./tests_tmp/km_dir_test/counts/partition_1/D1.kmer");
    size_t n = 0;
    while(kr.read<MK, MC>(kmer, count))
    {
      n++;
    }
    EXPECT_EQ(n, 46);
  }
  {
    km::Kmer<MK> kmer; kmer.set_k(31);
    uint32_t count;
    km::KmerReader<8192> kr("./tests_tmp/km_dir_test/counts/partition_2/D1.kmer");
    size_t n = 0;
    while(kr.read<MK, MC>(kmer, count))
    {
      n++;
    }
    EXPECT_EQ(n, 12);
  }
  {
    km::Kmer<MK> kmer; kmer.set_k(31);
    uint32_t count;
    km::KmerReader<8192> kr("./tests_tmp/km_dir_test/counts/partition_3/D1.kmer");
    size_t n = 0;
    while(kr.read<MK, MC>(kmer, count))
    {
      n++;
    }
    EXPECT_EQ(n, 43);
  }
  {
    km::Kmer<MK> kmer; kmer.set_k(31);
    uint32_t count;
    km::KmerReader<8192> kr("./tests_tmp/km_dir_test/counts/partition_0/D2.kmer");
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "AATATTATATCTACTACCATCATCATCACTA"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "AAGGAATATTATATCTACTACCATCATCATC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ATATCTTCCTCTCTTCGGGGGGGGGGGGGGG"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ATATTATATCTTCCTCTCTTCGGGGGGGGGG"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ATTATATCTTCCTCTCTTCGGGGGGGGGGGG"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ATGATGATGGTAGTAGATATAATATTCCTTC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "ATGATGGTAGTAGATATAATATTCCTTCCTC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "AGAGGAAGGAATATTATATCTACTACCATCA"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "AGTGATGATGATGGTAGTAGATATAATATTC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "AGGAATATTATATCTACTACCATCATCATCA"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "AGGAAGGAATATTATATCTACTACCATCATC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "CCCCCCCCCCCCCCGAAGAGAGGAAGATATA"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "CCCCCCCCCCCCCGAAGAGAGGAAGATATAA"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "CCCCCCCCCCCGAAGAGAGGAAGATATAATA"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "CCCCCCCCCGAAGAGAGGAAGATATAATATA"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "CCGCGTTTTTTTTTTTTTTTTTTTTCCCCCC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "TGATGATGGTAGTAGATATAATATTCCTTCC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "GCCGCGTTTTTTTTTTTTTTTTTTTTCCCCC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "GCGTTTTTTTTTTTTTTTTTTTTCCCCCCCC"); EXPECT_EQ(1, count);
    kr.read<MK, MC>(kmer, count);
    EXPECT_EQ(kmer.to_string(), "GTGATGATGATGGTAGTAGATATAATATTCC"); EXPECT_EQ(1, count);
    bool f = kr.read<MK, MC>(kmer, count); EXPECT_FALSE(f);
  }
  {
    km::Kmer<MK> kmer; kmer.set_k(31);
    uint32_t count;
    km::KmerReader<8192> kr("./tests_tmp/km_dir_test/counts/partition_1/D2.kmer");
    size_t n = 0;
    while(kr.read<MK, MC>(kmer, count))
    {
      n++;
    }
    EXPECT_EQ(n, 21);
  }
  {
    km::Kmer<MK> kmer; kmer.set_k(31);
    uint32_t count;
    km::KmerReader<8192> kr("./tests_tmp/km_dir_test/counts/partition_2/D2.kmer");
    size_t n = 0;
    while(kr.read<MK, MC>(kmer, count))
    {
      n++;
    }
    EXPECT_EQ(n, 58);
  }
  {
    km::Kmer<MK> kmer; kmer.set_k(31);
    uint32_t count;
    km::KmerReader<8192> kr("./tests_tmp/km_dir_test/counts/partition_3/D2.kmer");
    size_t n = 0;
    while(kr.read<MK, MC>(kmer, count))
    {
      n++;
    }
    EXPECT_EQ(n, 39);
  }
}

TEST(count_task, hash_count_task)
{
  km::KmDir::get().init(dir, "", false);
  Storage* config_storage = StorageFactory(STORAGE_FILE).load(km::KmDir::get().m_config_storage);
  LOCAL(config_storage);
  Configuration config = Configuration();
  config.load(config_storage->getGroup("gatb"));
  km::HashWindow hw("./data/hash.info");

  {
    km::sk_storage_t storage = std::make_shared<km::SuperKStorageReader>(km::KmDir::get().get_superk_path("D1"));
    km::parti_info_t pinfo = std::make_shared<PartiInfo<5>>(km::KmDir::get().get_superk_path("D1"));
    for (size_t p=0; p<4; p++)
    {
      std::string path = km::KmDir::get().get_count_part_path("D1", p, false, km::KM_FILE::HASH);
      km::HashCountTask<MK, MC, km::SuperKStorageReader> task(
        path, config, storage, pinfo, p, 0, hw.get_window_size_bits(), 31, 1, false, nullptr, false);
      task.exec();
    }
  }
  {
    km::sk_storage_t storage = std::make_shared<km::SuperKStorageReader>(km::KmDir::get().get_superk_path("D2"));
    km::parti_info_t pinfo = std::make_shared<PartiInfo<5>>(km::KmDir::get().get_superk_path("D2"));
    for (size_t p=0; p<4; p++)
    {
      std::string path = km::KmDir::get().get_count_part_path("D2", p, false, km::KM_FILE::HASH);
      km::HashCountTask<MK, MC, km::SuperKStorageReader> task(
        path, config, storage, pinfo, p, 0, hw.get_window_size_bits(), 31, 1, false, nullptr, false);
      task.exec();
    }
  }
  {
    uint64_t kmer;
    uint32_t count;
    km::HashReader<4294967296> kr("./tests_tmp/km_dir_test/counts/partition_0/D1.hash");
    kr.read(kmer, count); EXPECT_EQ(kmer, 155248); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 2023705); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 2567452); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 3271445); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 3722264); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 3868850); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 3981633); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 4503227); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 4962163); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 6435533); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 6862965); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 6978078); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 6979593); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 7059918); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 7083145); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 7093738); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 7725591); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 9582574); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 10836088); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 12171240); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 12224316); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 13314627); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 14513366); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 14877205); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 15672741); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 16842616); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 16978940); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 17200308); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 18924300); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 20011491); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 20323485); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 22465575); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 22637986); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 22862427); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 22918283); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 23401230); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 24946865); EXPECT_EQ(count, 1);
    bool f = kr.read(kmer, count); EXPECT_FALSE(f);
  }
  {
    uint64_t kmer;
    uint32_t count;
    km::HashReader<4294967296> kr("./tests_tmp/km_dir_test/counts/partition_1/D1.hash");
    size_t n = 0;
    while(kr.read(kmer, count))
    {
      n++;
    }
    EXPECT_EQ(n, 46);
  }
  {
    uint64_t kmer;
    uint32_t count;
    km::HashReader<4294967296> kr("./tests_tmp/km_dir_test/counts/partition_2/D1.hash");
    size_t n = 0;
    while(kr.read(kmer, count))
    {
      n++;
    }
    EXPECT_EQ(n, 12);
  }
  {
    uint64_t kmer;
    uint32_t count;
    km::HashReader<4294967296> kr("./tests_tmp/km_dir_test/counts/partition_3/D1.hash");
    size_t n = 0;
    while(kr.read(kmer, count))
    {
      n++;
    }
    EXPECT_EQ(n, 43);
  }
  {
    uint64_t kmer;
    uint32_t count;
    km::HashReader<4294967296> kr("./tests_tmp/km_dir_test/counts/partition_0/D2.hash");
    kr.read(kmer, count); EXPECT_EQ(kmer, 1303048); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 2821956); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 3573446); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 4954576); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 5940341); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 5964929); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 8761973); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 12178217); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 13532002); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 16524943); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 18299923); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 18309679); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 18709087); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 20543310); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 20906898); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 21688335); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 22116393); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 23796973); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 24160588); EXPECT_EQ(count, 1);
    kr.read(kmer, count); EXPECT_EQ(kmer, 24544513); EXPECT_EQ(count, 1);
    bool f = kr.read(kmer, count); EXPECT_FALSE(f);
  }
  {
    uint64_t kmer;
    uint32_t count;
    km::HashReader<4294967296> kr("./tests_tmp/km_dir_test/counts/partition_1/D2.hash");
    size_t n = 0;
    while(kr.read(kmer, count))
    {
      n++;
    }
    EXPECT_EQ(n, 21);
  }
  {
    uint64_t kmer;
    uint32_t count;
    km::HashReader<4294967296> kr("./tests_tmp/km_dir_test/counts/partition_2/D2.hash");
    size_t n = 0;
    while(kr.read(kmer, count))
    {
      n++;
    }
    EXPECT_EQ(n, 58);
  }
  {
    uint64_t kmer;
    uint32_t count;
    km::HashReader<4294967296> kr("./tests_tmp/km_dir_test/counts/partition_3/D2.hash");
    size_t n = 0;
    while(kr.read(kmer, count))
    {
      n++;
    }
    EXPECT_EQ(n, 39);
  }
}