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

#include <iostream>
#include <string>
#include <map>
#include <vector>

namespace km 
{

enum levels {
    DEBUG,
    INFO,
    WARN,
    ERROR
};

#ifdef _KM_LIB_INCLUDE_
extern const std::map<levels, std::string> level_to_string;
extern const std::map<std::string, levels> string_to_level;
#endif

#ifndef _KM_LIB_INCLUDE_
std::map<levels, std::string> level_to_string {
    {DEBUG, "DEBUG"},
    {INFO,  "INFO"},
    {WARN,  "WARN"},
    {ERROR, "ERROR"}
};

std::map<std::string, levels> string_to_level {
    {"DEBUG", DEBUG},
    {"INFO",  INFO},
    {"DEBUG", WARN},
    {"ERROR", ERROR}
};
#endif

struct log_config {
    bool show_labels = false;
    levels level = WARN;
};

extern log_config LOG_CONFIG;

class LOG
{
public:
    LOG (levels level, std::ostream&os = std::cerr)
        : level(level), cond(true), was_used(false), out(os) { init(); }

    LOG (levels level, bool conditional, std::ostream&os = std::cerr)
        : level(level), cond(conditional), was_used(false), out(os) { init(); }

    void init()
    {
        if (cond)
            if (level >= LOG_CONFIG.level)
            {
                if (LOG_CONFIG.show_labels)
                    operator << ("["+level_to_string.at(level)+"] - ");
                was_used = true;
            }
    }

    ~LOG ()
    {
        if (was_used)
        {
            if (level == DEBUG)
                out << std::endl;
            else 
                out << "\n";
            was_used = false;
        }
    }

    template<typename T>
    LOG &operator <<(const T& message)
    {
        if (cond)
            if (level >= LOG_CONFIG.level)
            {
                out << message;
                was_used = true;
            }
        return *this;
    }

    template<typename T>
    LOG &operator <<(const std::vector<T>& vector)
    {
        if (cond)
            if (level >= LOG_CONFIG.level)
            {
                for (auto &elem : vector)
                    out << elem << " ";
                was_used = true;
            }
        return *this;
    }

private:
    levels       level;
    bool         cond;
    bool         was_used;
    std::ostream &out;
};

} // namespace km 