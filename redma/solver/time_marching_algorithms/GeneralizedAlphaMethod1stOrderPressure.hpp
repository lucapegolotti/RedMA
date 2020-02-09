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

#ifndef GENERALIZEDALPHAMETHOD1STORDERPRESSURE_HPP
#define GENERALIZEDALPHAMETHOD1STORDERPRESSURE_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/time_marching_algorithms/aTimeMarchingAlgorithm.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/solver/array/BlockVector.hpp>
#include <redma/solver/system_solver/FunctionFunctor.hpp>
#include <redma/solver/system_solver/SystemSolver.hpp>
#include <redma/solver/array/Double.hpp>

#include <memory>

// for the details:
// https://arxiv.org/pdf/2002.01098.pdf
// In this modified version we use 1st order integration for pressure, as done
// in simvascular
namespace RedMA
{

template <class InVectorType, class InMatrixType>
class GeneralizedAlphaMethod1stOrderPressure : public aTimeMarchingAlgorithm<InVectorType, InMatrixType>
{
    typedef aFunctionProvider<InVectorType COMMA InMatrixType>  FunProvider;

public:

    GeneralizedAlphaMethod1stOrderPressure(const DataContainer& data);

    GeneralizedAlphaMethod1stOrderPressure(const DataContainer& data, SHP(FunProvider) funProvider);

    virtual void setup(const BlockVector<InVectorType>& zeroVector) override;

    virtual BlockVector<InVectorType> advance(const double& time, double& dt,
                                              int& status) override;

    virtual BlockVector<InVectorType> computeDerivative(const BlockVector<InVectorType>& solnp1,
                                                        double& dt) override;

    virtual void shiftSolutions(const BlockVector<InVectorType>& sol) override;

protected:
    BlockVector<InVectorType> computesolnp1(BlockVector<InVectorType> dersol,
                                            const double& dt);

    BlockVector<InVectorType> computesolnpalphaf(BlockVector<InVectorType> solnp1);

    BlockVector<InVectorType> computedersolnpalpham(BlockVector<InVectorType> dersol);


    BlockVector<InVectorType>              M_prevSolution;
    BlockVector<InVectorType>              M_prevDerivative;
    unsigned int                           M_order;
    double                                 M_alpham;
    double                                 M_alphaf;
    double                                 M_gamma;
    double                                 M_rhoinf;
    unsigned int                           M_nPrimalBlocks;
};

}

#include "GeneralizedAlphaMethod1stOrderPressure_imp.hpp"

#endif // GENERALIZEDALPHAMETHOD1STORDERPRESSURE_HPP
