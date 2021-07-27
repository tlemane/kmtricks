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
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

#include <kmtricks/io/fof.hpp>

namespace km {

void format_result_vector(const std::string& path,
                          std::ostream& stream,
                          std::vector<std::string>& query_idx,
                          Fof& fof);
void format_result_list(const std::string& path,
                        std::ostream& stream,
                        std::vector<std::string>& query_idx,
                        Fof& fof);

};