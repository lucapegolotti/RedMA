// Reduced Modeling of Arteries (RedMA)
// Copyright (C) 2019  Luca Pegolotti
//
// RedMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RedMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef PRINTLOG_HPP
#define PRINTLOG_HPP

#include <iostream>
#include <string>
#include <sstream>

namespace RedMA
{

enum Color {BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE};

extern void printlog(Color outColor, std::string text,
                     bool verbose = true);

template<typename T>
extern std::string to_string(const T& n);

extern void printlog(Color outColor, int num, bool verbose = true);

class CoutRedirecter
{
public:
    void redirect();

    std::string restore();

private:
    std::streambuf* M_prevBuf;
    std::ostringstream M_strCout;
};

}  // namespace RedMA

#endif  // PRINTLOG_HPP
