#include <gtest/gtest.h>
#include <filesystem>

namespace fs = std::filesystem;

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  fs::create_directory("./tests_tmp");
  int r = RUN_ALL_TESTS();
  //fs::remove_all("./tests_tmp/km_dir_test");
  return r;
}