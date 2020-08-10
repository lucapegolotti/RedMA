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
#include <redma/array/BlockVector.hpp>
#include <redma/solver/system_solver/SystemSolver.hpp>

#include <redma/problem/DataContainer.hpp>

#include <fstream>

namespace RedMA
{

class aTimeMarchingAlgorithm
{
    typedef aFunctionProvider         FunProvider;
public:

    aTimeMarchingAlgorithm(const DataContainer& datafile);

    aTimeMarchingAlgorithm(const DataContainer& datafile,
                           SHP(FunProvider) funProvider);

    virtual void setup(const BlockVector& zeroVector) = 0;

    virtual BlockVector advance(const double& time, double& dt,
                                int& status) = 0;

    // compute derivative of u at tn+1 given its value
    virtual BlockVector computeDerivative(const BlockVector& solnp1,
                                          double& dt) = 0;

    virtual void shiftSolutions(const BlockVector& sol) = 0;

    void dumpSolverStatistics(std::vector<SolverStatistics> statistics,
                              const double& t) const;

protected:
    void initializeStatisticsFile();

    DataContainer                                       M_data;
    SystemSolver                                        M_systemSolver;
    SHP(FunProvider)                                    M_funProvider;
    std::string                                         M_statisticsFile;
};

}

#endif // aTIMEMARCHINGALGORITHM_HPP
