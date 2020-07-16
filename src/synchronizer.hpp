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
#include <libgen.h>

#include "config.hpp"
#include <gatb/gatb_core.hpp>
#include <fmt/format.h>
#include <csignal>
#include <unistd.h>

class FBasedSync
{
public:
  FBasedSync();
  FBasedSync(
    vector<string>    banks,
    string            dir_sync,
    string            temp,
    size_t*           jobs,
    size_t            max_job,
    IteratorListener* progress,
    Env*              e,
    int               m
  );
  ~FBasedSync();
  void wait_ressources(size_t curr);
  void wait_end();
  void push(string s);
  static void write_end_signal(string& pref, string& synchro_dir);

private:
  vector<string>  _queue;
  vector<string>  _end_queue;
  vector<string>  _banks;
  size_t*         _jobs;
  size_t          _max_job;
  string          _dir_sync;
  string          _temp;
  IteratorListener* _progress;
  Env*              _e;
  int               _m;
};

FBasedSync::FBasedSync(
  vector<string>    banks,
  string            dir_sync,
  string            temp,
  size_t*           jobs,
  size_t            max_job,
  IteratorListener* progress,
  Env*              e,
  int               m
  ) : 
  _banks(banks),
  _dir_sync(dir_sync),
  _temp(temp),
  _max_job(max_job),
  _progress(progress),
  _jobs(jobs),
  _e(e),
  _m(m)
{
  _progress->init();
}

FBasedSync::~FBasedSync() {}

void FBasedSync::push(string s)
{
  _queue.push_back(s);
}

void FBasedSync::wait_ressources(size_t curr)
{
  if (*_jobs >= _max_job) {
    while (true)
    {
      bool job_available = false;
      for (size_t j=0; j<_queue.size(); j++) 
      {
        string pref = basename((char*)_queue[j].c_str());
        string end_signal = _dir_sync + fmt::format(_temp, pref);
        if (System::file().doesExist(end_signal))
        {
          _end_queue.push_back(_queue[j]);
          job_available = true;
          --*_jobs;
          _progress->inc(1);
        }
      }
      if (job_available)
      {
        for (size_t j=0; j<_end_queue.size(); j++)
        {
          _queue.erase(remove(_queue.begin(), _queue.end(), _end_queue[j]), _queue.end());
          if (_m)
          {
            string dpath = _e->STORE_KMERS + fmt::format(PART_DIR, _end_queue[j]);
            vector<string> dpaths = System::file().listdir(dpath);
            for (auto& path: dpaths)
              if (path.size()>2)
                System::file().remove(dpath + "/" + path);
          }

        }

        _end_queue.clear();
        break;
      }

      if (curr >= _banks.size()) break;
    }
  }
}

void FBasedSync::wait_end()
{
  while (*_jobs > 0)
  {
    bool finished_jobs = false;
    for (size_t j=0; j<_queue.size(); j++)
    {
      string pref = basename((char*)_queue[j].c_str());
      string end_signal = _dir_sync + fmt::format(_temp, pref);
      if (System::file().doesExist(end_signal))
      {
        _end_queue.push_back(_queue[j]);
        finished_jobs = true;
        --*_jobs;
        _progress->inc(1);
      }
    }

    if (finished_jobs)
    {
      for (size_t j=0; j<_end_queue.size(); j++)
      {
        _queue.erase(remove(_queue.begin(), _queue.end(), _end_queue[j]), _queue.end());
        if (_m)
        {
          string dpath = _e->STORE_KMERS + fmt::format(PART_DIR, _end_queue[j]);
          vector<string> dpaths = System::file().listdir(dpath);
          for (auto& path: dpaths)
            if (path.size()>2)
              System::file().remove(dpath + "/" + path);
        }
      }
      _end_queue.clear();
    }
  }
  _progress->finish();
  delete _progress;
}

void FBasedSync::write_end_signal(string& pref, string& synchro_dir)
{
  IFile* sync_file = System::file().newFile(synchro_dir+"/"+pref, "w");
  sync_file->flush();
  delete sync_file;
}

