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

#ifndef WOMERSLEYSOLUTION_H
#define WOMERSLEYSOLUTION_H

#include <AbstractFunctor.hpp>
#include "Womersley.hpp"

namespace RedMA
{

class WomersleySolution : public AbstractFunctor
{
    typedef std::function<double(double const&,
                                 double const&,
                                 double const&,
                                 double const&,
                                 unsigned int const& )>    Function;
public:
    ~WomersleySolution()
    {

    }

    double operator()(const double& t, const double& x, const double& y,
                      const double& z, const unsigned int& ic) const override
    {
        return LifeV::Womersley::uexact(t, x, y, z, ic);
    }

    double grad(unsigned int icoor, const double& t, const double& x, const double& y,
                const double& z, const unsigned int& ic) const override
    {
        return LifeV::Womersley::grad_u(icoor, t, x, y, z, ic);
    }

    Function exactFunction(const unsigned int& index) override
    {
        Function retFunction;
        if (index == 0)
            retFunction = LifeV::Womersley::uexact;
        else if (index == 1)
            retFunction = LifeV::Womersley::pexact;
        return retFunction;
    }
};

} //namespace RedMA

#endif  // WOMERSLEYSOLUTION_H_
