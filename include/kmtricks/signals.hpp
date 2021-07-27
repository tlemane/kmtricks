/*****************************************************************************
 *   kmdiff
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
#include <cxxabi.h>
#include <dlfcn.h>
#include <execinfo.h>

#include <csignal>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

#include <fmt/format.h>
#include <spdlog/spdlog.h>
#include <kmtricks/config.hpp>
#include <kmtricks/io/io_common.hpp>

namespace km
{
inline std::string signal_to_string(int signal)
{
  switch (signal)
  {
    case SIGABRT:
      return "SIGABRT";
      break;
    case SIGFPE:
      return "SIGFPE";
      break;
    case SIGILL:
      return "SIGILL";
      break;
    case SIGINT:
      return "SIGINT";
      break;
    case SIGSEGV:
      return "SIGSEGV";
      break;
    case SIGTERM:
      return "SIGTERM";
      break;
    default:
      return "?";
      break;
  }
}

class SignalHandler
{
 public:
  static SignalHandler& get()
  {
    static SignalHandler singleton;
    return singleton;
  }

  template <typename Callback = void (*)(int)>
  void init(Callback c = nullptr)
  {
    std::signal(SIGABRT, c ? c : default_callback);
    std::signal(SIGFPE, c ? c : default_callback);
    std::signal(SIGILL, c ? c : default_callback);
    std::signal(SIGINT, c ? c : default_callback);
    std::signal(SIGSEGV, c ? c : default_callback);
    std::signal(SIGTERM, c ? c : default_callback);
  }

  template <typename Callback = void (*)(int)>
  void set(int signal, Callback c)
  {
    std::signal(signal, c);
  }

  static void default_callback(int signal)
  {
    const std::string str_signal = strsignal(signal);
    std::stringstream ss;
    int size;
    void* stack[256];
    int max_size = sizeof(stack) / sizeof(stack[0]);
    char buffer[1024];

    size = backtrace(stack, max_size);
    char** symbols = backtrace_symbols(stack, size);

    std::string backtrace_path = fmt::format("./{}_backtrace.log", PROJECT_NAME);
    ss << "\nBacktrace:\n";
    for (int i = 1; i < size; i++)
    {
      Dl_info info;
      if (dladdr(stack[i], &info))
      {
        char* name = NULL;
        int status;
        name = abi::__cxa_demangle(info.dli_sname, NULL, 0, &status);

        ss << i << " ";
        ss << "0x" << std::setfill('0') << std::setw(sizeof(void*) * 2) << std::hex;
        ss << reinterpret_cast<uint64_t>(stack[i]) << " " << std::dec;
        snprintf(buffer, sizeof(buffer), "%s", status == 0 ? name : info.dli_sname);
        ss << buffer;
        ss << " + " << static_cast<char*>(stack[i]) - static_cast<char*>(info.dli_saddr) << "\n";
        free(name);
      }
      else
      {
        ss << i << " " << stack[i];
      }
    }
    free(symbols);

    if (size == max_size) ss << "[truncated]\n";

    if (signal != SIGINT)
    {
      std::ofstream signal_file(backtrace_path, std::ios::out);
      signal_file << ss.str();
      signal_file.close();
    }
    std::stringstream msg;
    msg << fmt::format(
        "Killed after receive {}:{}({}) signal.", str_signal, signal_to_string(signal), signal);
    if (signal != SIGINT)
    {
      msg << fmt::format(" Demangled backtrace dumped at {}. ", backtrace_path);
      const std::string issue =
          "If the problem persists, please open an issue with the return of '{} infos' "
          "and the content of {}";
      msg << fmt::format(issue, PROJECT_NAME, backtrace_path);
    }
    spdlog::error(msg.str());

    exit(signal);
  }

 private:
  SignalHandler() {}
};

};  // namespace kmdiff