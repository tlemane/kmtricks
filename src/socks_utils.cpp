#include <kmtricks/utils.hpp>
#include <kmtricks/socks-interface/socks_utils.hpp>

namespace km {

void format_result_vector(const std::string& path,
                          std::ostream& stream,
                          std::vector<std::string>& query_idx,
                          Fof& fof)
{
  std::ifstream in(path, std::ios::in); check_fstream_good(path, in);
  std::vector<char> res(fof.size(), '0');
  std::string name = "";
  size_t qid = 0;
  for (std::string line; std::getline(in, line);)
  {
    if (line[0] == '*')
    {
      if (!name.empty())
      {
        stream << name << ": ";
        for (auto& c : res)
          stream << c << ' ';
        stream << "\n";
        std::fill(res.begin(), res.end(), '0');
        name = query_idx[qid]; qid++;
      }
      else
      {
        name = query_idx[qid]; qid++;
      }
    }
    else
    {
      res[fof.get_i(line)] = '1';
    }
  }
  stream << name << ": ";
  for (auto& c : res)
    stream << c << ' ';
  stream << "\n";
}

void format_result_list(const std::string& path,
                        std::ostream& stream,
                        std::vector<std::string>& query_idx,
                        Fof& fof)
{
  std::ifstream in(path, std::ios::in); check_fstream_good(path, in);
  std::vector<std::string> res;
  std::string name = "";
  size_t qid = 0;
  for (std::string line; std::getline(in, line);)
  {
    if (line[0] == '*')
    {
      if (!name.empty())
      {
        stream << name << ": ";
        for (auto& s : res)
          stream << s << ' ';
        stream << "\n";
        res.clear();
        name = query_idx[qid]; qid++;
      }
      else
      {
        name = query_idx[qid]; qid++;
      }
    }
    else
    {
      res.push_back(line);
    }
  }
  stream << name << ": ";
  for (auto& s : res)
    stream << s << ' ';
  stream << "\n";
}

};