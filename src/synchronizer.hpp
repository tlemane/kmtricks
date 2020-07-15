#ifndef _SYNCHRONIZER_HPP_
#define _SYNCHRONIZER_HPP_

#include <libgen.h>

#include "config.hpp"
#include <gatb/gatb_core.hpp>
#include <fmt/format.h>
#include <csignal>
#include <unistd.h>

typedef string    cmd_t;
typedef string    kbf_t;
typedef string    sync_t;
typedef uint32_t  rmem_t;

typedef map<cmd_t, tuple<sync_t, rmem_t>>Commands;
typedef map<pid_t, tuple<kbf_t, sync_t>> pids_t;
typedef tuple<cmd_t, sync_t, rmem_t> job_t;

//class FBasedSync2
//{
//public:
//  FBasedSync2(
//    Commands &cmds,
//    Env      *env,
//    size_t   max_cores,
//    size_t   max_mem,
//    IteratorListener* progress
//    );
//  void execute();
//
//private:
//  void check();
//  tuple<cmd_t, sync_t, rmem_t> find();
//  void wait_ressource();
//  void execj(job_t &j);
//
//private:
//  Commands _cmds;
//  Env      *_e;
//  size_t   _max_cores;
//  size_t   _max_mem;
//  size_t   _av_cores;
//  size_t   _av_mem;
//  size_t    _nb_cmds;
//  pids_t   _current;
//  pids_t   _finish;
//  size_t   av_jobs;
//  IteratorListener* _progress;
//};
//
//tuple<cmd_t, sync_t, rmem_t> FBasedSync2::find()
//{
//  job_t job = make_tuple("", "", 0);
//  for ( auto &cmd: _cmds )
//  {
//    rmem_t rmem = get<1>(cmd.second);
//    if ( rmem < _av_mem )
//    {
//      sync_t tmp_sync = get<0>(cmd.second);
//      _cmds.erase(cmd.first);
//      job = make_tuple(cmd.first, tmp_sync, rmem);
//      break;
//    }
//  }
//  if (get<2>(job))
//    _cmds.erase(get<0>(job));
//  return job;
//}
//
//FBasedSync2::FBasedSync2(Commands &cmds, Env *env, size_t max_cores, size_t max_mem, IteratorListener* progress)
//  : _cmds(cmds), _e(env), _max_cores(max_cores), _max_mem(max_mem), _progress(progress)
//{
//  _av_cores = _max_cores;
//  _av_mem   = _max_mem;
//  _nb_cmds  = _cmds.size();
//}
//
//void FBasedSync2::wait_ressource()
//{
//
//}
//
//void FBasedSync2::execj(job_t &j)
//{
//  av_jobs--;
//  pid_t fpid = fork();
//  if (fpid == 0)
//  {
//    pid_t cpid = getpid();
//    _current[cpid] = make_tuple("", get<1>(j));
//    execv("A", get<0>(j).c_str());
//  }
//  else
//  {
//    _av_mem += get<2>(j);
//  }
//}
//void FBasedSync2::execute()
//{
//  size_t max_concurrent = min(_max_mem/mem_p_j, _max_cores);
//  av_jobs = max_concurrent;
//  size_t _nb_run = 0;
//  while (_nb_cmds != _nb_run)
//  {
//    if (!av_jobs)
//      wait_ressource();
//    else
//      job_t new_job = find();
//  }
//  if (av_jobs)
//  {
//    job_t new_j = find();
//  }
//  for (auto& cmd: _cmds)
//  {
//    cmd_t   curr_command = cmd.first;
//    sync_t  curr_fsync = get<0>(cmd.second);
//    rmem_t  rmem       = get<1>(cmd.second);
//    if (av_jobs > 0)
//    {
//      pid_t fpid = fork();
//      if (fpid == 0)
//      {
//        pid_t cpid = getpid();
//        _current[cpid] = make_tuple("kbf_couter", curr_fsync);
//      }
//      else {} //
//    }
//    else
//    {
//
//    }
//  }
//}
//
//void FBasedSync2::check()
//{
//  for (auto& child: _current)
//  {
//    pid_t cpid = child.first;
//    kbf_t  kbf_bin = get<0>(child.second);
//    sync_t fsync = get<1>(child.second);
//
//    if (0 != kill(cpid, 0))
//      if (!System::file().doesExist(fsync))
//      {
//        cout << "Process " << fsync << " stopped before sending its end signal" << endl;
//        std::system(fmt::format("killall {}", kbf_bin).c_str());
//        cout << "Execution halted" << endl;
//        exit(EXIT_FAILURE);
//      }
//  }
//}

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

#endif // _SYNCHRONIZER_HPP_