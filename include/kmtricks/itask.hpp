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
#include <functional>
#include <cstdint>
#include <memory>

namespace km {

class ITask
{
public:
  ITask (uint32_t level, bool clear = false) : m_priority_level(level), m_clear(clear) {}

  virtual ~ITask() = default;

  virtual void preprocess() = 0;
  virtual void postprocess() = 0;
  virtual void exec() = 0;

  virtual void set_level(uint32_t level) { m_priority_level = level; }

  bool operator==(const ITask& task) const
  {
    return m_priority_level == task.m_priority_level;
  }

  bool operator>(const ITask& task) const
  {
    return m_priority_level > task.m_priority_level;
  }

  bool operator<(const ITask& task) const
  {
    return m_priority_level < task.m_priority_level;
  }

  void set_callback(std::function<void()> callback)
  {
    m_callback = callback;
  }

  void exec_callback()
  {
    if (m_callback && !m_ce)
    {
      m_callback();
      m_ce = true;
    }
  }

  bool finish() const
  {
    return m_finish;
  }

  bool running() const
  {
    return m_running;
  }

  bool in_queue() const
  {
    return m_in_queue;
  }

  void in() {m_in_queue = true;}
  void out() {m_in_queue = false;}

protected:
  uint32_t m_priority_level;
  bool m_ready {false};
  bool m_finish {false};
  bool m_ce {false};
  bool m_clear {false};
  bool m_running {false};
  bool m_in_queue {false};
  std::function<void()> m_callback {nullptr};
};

using task_t = std::shared_ptr<ITask>;


inline bool operator==(const task_t& lhs, const task_t& rhs)
{
  return (*lhs.get()) == (*rhs.get());
}

inline bool operator>(const task_t& lhs, const task_t& rhs)
{
  return (*lhs.get()) > (*rhs.get());
}

inline bool operator<(const task_t& lhs, const task_t& rhs)
{
  return (*lhs.get()) < (*rhs.get());
}

};