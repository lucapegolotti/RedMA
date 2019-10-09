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

#ifndef ABSTRACTFUNCTOR_H
#define ABSTRACTFUNCTOR_H

namespace RedMA
{

class AbstractFunctor
{
    typedef std::function<double(double const&,
                                 double const&,
                                 double const&,
                                 double const&,
                                 unsigned int const& )>    Function;

public:
    virtual ~AbstractFunctor(){};

    virtual double operator()(const double& t, const double& x, const double& y,
                              const double& z, const unsigned int& ic) const = 0;

    virtual double grad(unsigned int icoor, const double& t, const double& x, const double& y,
                        const double& z, const unsigned int& ic) const = 0;

    virtual Function exactFunction(const unsigned int& index) = 0;

    virtual Function exactFunctionDt(const unsigned int& index) = 0;
};
}
#endif
