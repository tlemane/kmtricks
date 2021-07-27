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
#include <chrono>
#include <functional>
#include <string>
#include <string>
#include <sstream>
#include <iomanip>

namespace km {

class Timer
{
  using time_point_t = std::chrono::time_point<std::chrono::steady_clock>;
  using days = std::chrono::duration<int, std::ratio<86400>>;

 public:
  Timer() {start();}

  template <typename Unit>
  auto elapsed()
  {
    if (m_running) end();
    return std::chrono::duration_cast<Unit>(m_end_time - m_start_time);
  }

  std::string formatted()
  {
    if (m_running) end();
    std::chrono::seconds seconds(std::chrono::duration_cast<std::chrono::seconds>(m_end_time - m_start_time));
    auto d = std::chrono::duration_cast<days>(seconds); seconds -= d;
    auto h = std::chrono::duration_cast<std::chrono::hours>(seconds); seconds -= h;
    auto m = std::chrono::duration_cast<std::chrono::minutes>(seconds); seconds -= m;
    auto s = std::chrono::duration_cast<std::chrono::seconds>(seconds);

    std::stringstream ss; ss.fill('0');
    if (d.count())
      ss << std::setw(2) << d.count() << "d";
    if (h.count())
      ss << std::setw(2) << h.count() << "h";
    if (m.count())
      ss << std::setw(2) << m.count() << "m";
    ss << std::setw(2) << s.count() << "s";
    return ss.str();
  }

  template <typename Unit>
  static auto time_it(std::function<void()> func)
  {
    auto start = std::chrono::steady_clock::now();
    func();
    auto end = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<Unit>(end - start).count();
  }

  void reset()
  {
    m_running = false;
    m_start_time = time_point_t{};
    m_end_time = time_point_t{};
    start();
  }

 private:
  void start()
  {
    m_running = true;
    m_start_time = std::chrono::steady_clock::now();
  }
  void end()
  {
    m_running = false;
    m_end_time = std::chrono::steady_clock::now();
  }

 private:
  time_point_t m_start_time{};
  time_point_t m_end_time{};
  bool m_running{false};
};

};

