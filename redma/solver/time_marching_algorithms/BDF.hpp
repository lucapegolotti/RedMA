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

#ifndef BDF_HPP
#define BDF_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/time_marching_algorithms/aTimeMarchingAlgorithm.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/solver/array/BlockVector.hpp>
#include <redma/solver/system_solver/FunctionFunctor.hpp>
#include <redma/solver/system_solver/SystemSolver.hpp>
#include <redma/solver/array/Double.hpp>

#include <memory>

namespace RedMA
{

template <class InVectorType, class InMatrixType>
class BDF : public aTimeMarchingAlgorithm<InVectorType, InMatrixType>
{
    typedef aFunctionProvider<InVectorType COMMA InMatrixType>  FunProvider;

public:

    BDF(const DataContainer& data);

    BDF(const DataContainer& data, SHP(FunProvider) funProvider);

    virtual void setup(const BlockVector<InVectorType>& zeroVector) override;

    virtual BlockVector<InVectorType> advance(const double& time, double& dt,
                                              int& status) override;

    virtual void shiftSolutions(const BlockVector<InVectorType>& sol) override;

    virtual BlockVector<InVectorType> computeDerivative(const BlockVector<InVectorType>& solnp1,
                                                        double& dt) override;

    BlockVector<InVectorType> computeExtrapolatedSolution();

protected:
    std::vector<BlockVector<InVectorType>>   M_prevSolutions;
    std::vector<double>                      M_coefficients;
    unsigned int                             M_order;
    double                                   M_rhsCoeff;
    bool                                     M_useExtrapolation;
};

}

#include "BDF_imp.hpp"

#endif // BDF_HPP
