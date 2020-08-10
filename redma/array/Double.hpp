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

#ifndef DOUBLE_HPP
#define DOUBLE_HPP

#include <iostream>

namespace RedMA
{

// this clas is only to be used as template argument of block structures
class Double
{
public:
    Double();

    Double operator+(const Double& other);

    Double& operator+=(const Double& other);

    Double& operator-=(const Double& other);

    Double& operator*=(const double& coeff);

    Double operator*(const Double& other);

    Double& operator=(const double& other);

    void hardCopy(const Double& other);

    void softCopy(const Double& other);

    double norm2() const;

    double& data();

    double data() const;

    void dump(std::string filename) const {}

    Double block(unsigned int i, unsigned int j) {return *this;}

private:
    double  M_double;
};

}

#endif // VECTOREP_HPP
