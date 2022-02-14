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
 *
 *  @file kmtricks's exceptions
 *****************************************************************************/

#pragma once

#include <stdexcept>
#include <string>

namespace km
{
class km_exception : public std::exception
{
 private:
  std::string name{"Base error"};
  std::string msg{"Base error msg should never be printed"};

 public:
  km_exception(const std::string& name, const std::string& msg) : name(name), msg(msg) {}
  std::string get_name() const { return name; }
  std::string get_msg() const { return msg; }
};

#define km_EXCEPTION(name)                                     \
  class name : public km_exception                             \
  {                                                            \
   public:                                                     \
    name(const std::string& msg) : km_exception(#name, msg) {} \
  }

km_EXCEPTION(IOError);
km_EXCEPTION(IDError);
km_EXCEPTION(InputError);
km_EXCEPTION(FileNotFoundError);
km_EXCEPTION(PipelineError);
km_EXCEPTION(ConfigError);
km_EXCEPTION(KSizeError);
km_EXCEPTION(PluginError);

};
