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

// std
#include <functional>
#include <future>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <thread>
#include <memory>
#include <vector>

#include <kmtricks/itask.hpp>

namespace km
{
class TaskPool
{
  using size_type = std::result_of<decltype (&std::thread::hardware_concurrency)()>::type;

 public:
  TaskPool(size_type threads)
  {
    if (threads < m_n) m_n = threads;
    for (size_t i = 0; i < m_n; i++)
    {
      m_pool.push_back(std::thread(&TaskPool::worker, this, i));
    }
  }

  ~TaskPool()
  {
    {
      std::unique_lock<std::mutex> lock(m_queue_mutex);
      m_stop = true;
    }
    m_condition.notify_all();
    for (std::thread& t : m_pool)
      if (t.joinable()) t.join();
  }

  TaskPool() = delete;
  TaskPool(const TaskPool&) = delete;
  TaskPool& operator=(const TaskPool&) = delete;
  TaskPool(TaskPool&&) = delete;
  TaskPool& operator=(TaskPool&&) = delete;


  void join_all()
  {
    {
      std::unique_lock<std::mutex> lock(m_queue_mutex);
      m_stop = true;
    }
    m_condition.notify_all();
    for (std::thread& t : m_pool) t.join();
  }

  void join(int i)
  {
    if (m_pool[i].joinable()) m_pool[i].join();
  }

  void add_task(task_t task)
  {
    {
      std::unique_lock<std::mutex> lock(m_queue_mutex);
      task->in();
      m_queue.push(task);
    }
    m_condition.notify_one();
  }

 private:
  void worker(int i)
  {
    while (true)
    {
      task_t task;
      {
        std::unique_lock<std::mutex> lock(this->m_queue_mutex);
        this->m_condition.wait(lock, [this] { return this->m_stop || !this->m_queue.empty(); });
        if (this->m_stop && this->m_queue.empty()) return;
        task = this->m_queue.top();
        this->m_queue.pop();
      }
      task->preprocess();
      task->exec();
      task->postprocess();
      task->out();
    }
  }

 private:
  size_type m_n{std::thread::hardware_concurrency()};
  std::vector<std::thread> m_pool;
  std::priority_queue<task_t> m_queue;
  std::mutex m_queue_mutex;
  std::condition_variable m_condition;
  bool m_stop{false};
};

};
