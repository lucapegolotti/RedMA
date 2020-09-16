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

#ifndef GENERALIZEDALPHAMETHOD_HPP
#define GENERALIZEDALPHAMETHOD_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/time_marching_algorithms/aTimeMarchingAlgorithm.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/array/BlockVector.hpp>
#include <redma/solver/system_solver/FunctionFunctor.hpp>
#include <redma/solver/system_solver/SystemSolver.hpp>
#include <redma/array/Double.hpp>

#include <memory>

// for the details:
// https://arxiv.org/pdf/2002.01098.pdf

namespace RedMA
{

class GeneralizedAlphaMethod : public aTimeMarchingAlgorithm
{
    typedef aFunctionProvider  FunProvider;

public:

    GeneralizedAlphaMethod(const DataContainer& data);

    GeneralizedAlphaMethod(const DataContainer& data, SHP(FunProvider) funProvider);

    virtual void setup(const SHP(aVector)& zeroVector) override;

    virtual SHP(aVector) advance(const double& time, double& dt, int& status) override;

    virtual SHP(aVector) computeDerivative(const SHP(aVector)& solnp1, double& dt) override;

    virtual void shiftSolutions(const SHP(aVector)& sol) override;

protected:
    SHP(aVector) computesolnp1(SHP(aVector) dersol, const double& dt);

    SHP(aVector) computesolnpalphaf(SHP(aVector) solnp1);

    SHP(aVector) computedersolnpalpham(SHP(aVector) dersol);


    SHP(aVector)                            M_prevSolution;
    SHP(aVector)                            M_prevDerivative;
    unsigned int                           M_order;
    double                                 M_alpham;
    double                                 M_alphaf;
    double                                 M_gamma;
    double                                 M_rhoinf;
};

}

#endif // GENERALIZEDALPHAMETHOD_HPP
