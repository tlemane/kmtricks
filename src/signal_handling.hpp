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
#include <csignal>
#include <cstring>
#include <string>
#include <cstdio>
#include <iostream>
#include <unistd.h>
#include <execinfo.h>
#include <libgen.h>
#include "config.hpp"

#if __APPLE__
#include <mach-o/dyld.h>
#endif

using namespace std;

#define INIT_SIGN \
    signal(SIGABRT, SignalHandler::callback); \
    signal(SIGFPE, SignalHandler::callback);  \
    signal(SIGILL, SignalHandler::callback);  \
    signal(SIGINT, SignalHandler::callback);  \
    signal(SIGSEGV, SignalHandler::callback); \
    signal(SIGTERM, SignalHandler::callback) \

vector<string> split_str(string s, string t)
{
  vector<string> r;
  while(!s.empty())
  {
    int i = s.find(t);
    if (i != string::npos)
    {
      r.push_back(s.substr(0, i));
      s = s.substr(i+t.size());
      if (s.empty()) r.push_back(s);
    }
    else
    {
      r.push_back(s);
      s = "";
    }
  }
  return r;
}

class SignalHandler
{
public:
  explicit SignalHandler(string exec_name = "", int file_id = 0, int part_id = 0)
    : _exec_name(exec_name), _fid(file_id), _pid(part_id)
  {

  };
  ~SignalHandler() = default;

  bool check_dir()
  {
    if (System::file().doesExistDirectory("./km_backtrace"))
    {
      auto paths = System::file().listdir("./km_backtrace");
      if ( paths.size() > 2 )
      {
        for ( auto &p: paths )
        {
          if ( p.size() > 2 )
          {
            if (p.find("-") != string::npos)
            {
              vector<string> info = split_str(p, "-");
              _exec_name = info[0];
              _type = info[1];
            }
          }
        }
        kill_all();
        return true;
      }
    }
    return false;
  };

  void kill_all()
  {
    cerr << endl;
    cerr << fmt::format(ERROR_MSG,
      _type,
      _exec_name,
      CONTACT,
      BACKTRACE);
    std::system(KILLALL.c_str());
    exit(EXIT_FAILURE);
  }

  static string signal_to_string(int sig)
  {
    switch ( sig )
    {
      case SIGABRT: return "SIGARBT"; break;
      case SIGFPE : return "SIGFPE"; break;
      case SIGILL : return "SIGILL"; break;
      case SIGINT : return "SIGINT"; break;
      case SIGSEGV : return "SIGSEV"; break;
      case SIGTERM: return "SIGTERM"; break;
      default: return "UNKNOW"; break;
    }
  };

  static void callback(int sig)
  {
    string type = signal_to_string(sig);
    char buffer[1024];

#if __APPLE__
    uint32_t size=1024;
    _NSGetExecutablePath(buffer, &size);
#else
    readlink("/proc/self/exe", buffer, 1024);
#endif
    string bin = basename(buffer);
    System::file().mkdir("./km_backtrace", -1);
    ofstream fout(fmt::format(RUN_INFOS, bin, type), ios::out);
    fout.flush();
    fout.close();
    FILE *back = fopen(BACKTRACE.c_str(), "w");
    void *array[256];
    size_t _size;
    _size = backtrace(array, 256);
    backtrace_symbols_fd(array, _size, fileno(back));
  }

private:
  string _exec_name;
  string _type;
  int _fid;
  int _pid;
};

