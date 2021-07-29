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
 *****************************************************************************/

#pragma once
#include <indicators/multi_progress.hpp>
#include <indicators/dynamic_progress.hpp>
#include <indicators/progress_bar.hpp>

using namespace indicators;

namespace km {

inline ProgressBar* get_progress_bar(const std::string& name, size_t size,
                                     size_t width, Color color, bool time)
{
  return new ProgressBar{option::BarWidth{width},
                     option::Start{"["},
                     option::Fill{"="},
                     option::Lead{">"},
                     option::Remainder{" "},
                     option::End{"]"},
                     option::ForegroundColor{color},
                     option::ShowElapsedTime{true},
                     option::ShowRemainingTime{false},
                     option::PrefixText{name},
                     option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},
                     option::MaxProgress{size}
  };
}

};