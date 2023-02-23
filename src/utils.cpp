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

#include <kmtricks/utils.hpp>

namespace km {

VerbosityLevel str_to_verbosity_level(const std::string& str_level)
{
  if (str_level == "debug")
    return VerbosityLevel::DEBUG;
  else if (str_level == "info")
    return VerbosityLevel::INFO;
  else if (str_level == "warning")
    return VerbosityLevel::WARNING;
  else if (str_level == "error")
    return VerbosityLevel::ERROR;
  else
    return VerbosityLevel::WARNING;
}

void set_verbosity_level(const std::string& level)
{
  switch (str_to_verbosity_level(level))
  {
    case VerbosityLevel::DEBUG:
      spdlog::set_level(spdlog::level::debug);
      break;
    case VerbosityLevel::INFO:
      spdlog::set_level(spdlog::level::info);
      break;
    case VerbosityLevel::WARNING:
      spdlog::set_level(spdlog::level::warn);
      break;
    case VerbosityLevel::ERROR:
      spdlog::set_level(spdlog::level::err);
      break;
  }
}

};
