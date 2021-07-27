#include <iostream>
#include <filesystem>
int main()
{
  std::cout << std::filesystem::temp_directory_path();
  return 0;
}