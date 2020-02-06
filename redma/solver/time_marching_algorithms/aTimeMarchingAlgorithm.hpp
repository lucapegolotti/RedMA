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

#ifndef aTIMEMARCHINGALGORITHM_HPP
#define aTIMEMARCHINGALGORITHM_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/time_marching_algorithms/aFunctionProvider.hpp>
#include <redma/solver/array/BlockVector.hpp>
#include <redma/solver/system_solver/SystemSolver.hpp>

#include <redma/solver/problem/DataContainer.hpp>

#include <fstream>

namespace RedMA
{

template <class InVectorType, class InMatrixType>
class aTimeMarchingAlgorithm
{
    typedef aFunctionProvider<InVectorType COMMA InMatrixType>  FunProvider;
public:

    aTimeMarchingAlgorithm(const DataContainer& datafile);

    aTimeMarchingAlgorithm(const DataContainer& datafile,
                           SHP(FunProvider) funProvider);

    virtual void setup(const BlockVector<InVectorType>& zeroVector) = 0;

    virtual BlockVector<InVectorType> advance(const double& time, double& dt,
                                              int& status) = 0;

    // compute derivative of u at tn+1 given its value
    virtual BlockVector<InVectorType> computeDerivative(const BlockVector<InVectorType>& solnp1,
                                                        double& dt) = 0;

    // this must be implemented by multistep methods (e.g. bdf)
    virtual void shiftSolutions(const BlockVector<InVectorType>& sol) {}

    void dumpSolverStatistics(std::vector<SolverStatistics> statistics,
                              const double& t) const;

protected:
    void initializeStatisticsFile();

    DataContainer                                       M_data;
    SystemSolver<InVectorType, InMatrixType>            M_systemSolver;
    SHP(FunProvider)                                    M_funProvider;
    std::string                                         M_statisticsFile;
};

}

#include "aTimeMarchingAlgorithm_imp.hpp"

#endif // aTIMEMARCHINGALGORITHM_HPP
