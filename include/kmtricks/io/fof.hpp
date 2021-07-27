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
#include <regex>
#include <fstream>
#include <tuple>
#include <set>
#include <unordered_map>
#include <filesystem>

#include <bcli/bcli.hpp>

#include <kmtricks/utils.hpp>

namespace fs = std::filesystem;

namespace km {

class Fof
{
public:
  inline static std::regex pattern
    {R"((^[A-Za-z0-9_-]+)[\s]*:[\s]*([.A-Za-z0-9\/_\-; ]+)([\s]*![\s]*)?([0-9]+$)?)"};

  inline static std::regex invalid
    {R"(([<>{},[\]]))"};

  using data_t = std::vector<std::tuple<std::string, std::vector<std::string>, uint32_t>>;

public:
  Fof () {}

  Fof(const std::string& path)
    : m_path(path)
  {
    parse();
  }

  std::string get_all()
  {
    std::vector<std::string> p;
    for (auto [_, paths, __] : m_data)
    {
      std::copy(paths.begin(), paths.end(), std::back_inserter(p));
    }
    if (p.size() == 1)
    {
      return p[0];
    }
    std::ostringstream ss;
    std::copy(p.begin(), p.end(), std::ostream_iterator<std::string>(ss, ","));
    return ss.str();
  }

  std::string get_id(uint32_t i)
  {
    return std::get<0>(m_data[i]);
  }

  uint32_t get_i(const std::string& id)
  {
    if (!m_map.count(id))
      throw IDError(fmt::format("Unknown id: {}", id));
    return m_map.at(id);
  }

  std::string get_files(const std::string& id)
  {
    if (!m_map.count(id))
      throw IDError(fmt::format("Unknown id: {}", id));
    return bc::utils::join(std::get<1>(m_data[m_map.at(id)]), ",");
  }

  void copy(const std::string& path)
  {
    fs::copy_file(m_path, path);
  }

  size_t size() const
  {
    return m_data.size();
  }

  size_t total() const
  {
    size_t count = 0;
    for (auto& id : m_data)
      count += std::get<1>(id).size();
    return count;
  }

  auto begin() { return m_data.begin(); }
  auto end() { return m_data.end(); }
  auto begin() const { return m_data.cbegin(); }
  auto end() const { return m_data.cend(); }

private:
  void parse()
  {
    std::ifstream in(m_path, std::ios::in); check_fstream_good(m_path, in);
    uint32_t idx = 0;
    for (std::string line; std::getline(in, line);)
    {
      if (bc::utils::trim(line).empty())
        continue;
      std::smatch g, i;
      bool found = std::regex_search(line, g, Fof::pattern);
      bool inv = std::regex_search(line, i, Fof::invalid);
      if (!found || inv)
        throw IOError("Invalid fof format.");
      if (m_id.count(g[1].str()))
        throw IOError(fmt::format("{} -> sample identifiers must be unique.", g[1].str()));
      m_id.insert(g[1].str());
      std::vector<std::string> paths = bc::utils::split(g[2].str(), ';');
      std::transform(paths.begin(), paths.end(), paths.begin(), [](const std::string& s){
        return bc::utils::trim(s);
      });

      m_data.emplace_back(g[1].str(), paths,
                          g[4].str().empty() ? 0 : bc::utils::lexical_cast<uint32_t>(g[4].str()));
      m_map.insert({std::get<0>(m_data[idx]), idx});
      idx++;
    }
  }

private:
  std::string m_path;
  data_t m_data;
  std::set<std::string> m_id;
  std::unordered_map<std::string, uint32_t> m_map;
};

};