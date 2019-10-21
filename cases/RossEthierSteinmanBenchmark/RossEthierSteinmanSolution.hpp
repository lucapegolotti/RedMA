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

#ifndef ROSSETHIERSTEINMANSOLUTION_H
#define ROSSETHIERSTEINMANSOLUTION_H

#include <AbstractFunctor.hpp>
#include "RossEthierSteinmanDec.hpp"

namespace RedMA
{

class RossEthierSteinmanSolution : public AbstractFunctor
{
    typedef std::function<double(double const&,
                                 double const&,
                                 double const&,
                                 double const&,
                                 unsigned int const& )>    Function;

     typedef std::function<double(double const&,
                                  double const&,
                                  double const&,
                                  double const&,
                                  unsigned int const&,
                                  LifeV::VectorSmall<3> const&)>    FunctionNeumann;
public:
    ~RossEthierSteinmanSolution()
    {

    }

    double operator()(const double& t, const double& x, const double& y,
                      const double& z, const unsigned int& ic) const override
    {
        return LifeV::RossEthierSteinmanDec::uexact(t, x, y, z, ic);
    }

    double grad(unsigned int icoor, const double& t, const double& x, const double& y,
                const double& z, const unsigned int& ic) const override
    {
        return LifeV::RossEthierSteinmanDec::grad_u(icoor, t, x, y, z, ic);
    }

    Function exactFunction(const unsigned int& index) override
    {
        Function retFunction;
        if (index == 0)
            retFunction = LifeV::RossEthierSteinmanDec::uexact;
        else if (index == 1)
            retFunction = LifeV::RossEthierSteinmanDec::pexact;
        return retFunction;
    }

    Function exactFunctionDt(const unsigned int& index) override
    {
        Function retFunction;
        if (index == 0)
            retFunction = LifeV::RossEthierSteinmanDec::uexact_dt;
        else if (index == 1)
            retFunction = LifeV::RossEthierSteinmanDec::pexact_dt;
        return retFunction;
    }

    FunctionNeumann exactNeumann() override
    {
        return LifeV::RossEthierSteinmanDec::fNeumann;
    }

    FunctionNeumann exactNeumannDt() override
    {
        return LifeV::RossEthierSteinmanDec::fNeumann_dt;
    }
};

} //namespace RedMA

#endif  // ROSSETHIERSTEINMANSOLUTION_H
