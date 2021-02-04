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
                           shp<FunProvider> funProvider);

    aTimeMarchingAlgorithm(const DataContainer& datafile,
                           const shp<aVector>& zeroVector);

    virtual void setup(const shp<aVector>& zeroVector) = 0;

    virtual shp<aVector> advance(const double& time, double& dt,
                                int& status) = 0;

    virtual shp<aVector> simpleAdvance(const double &dt, const shp<BlockVector> &sol) = 0;

    // compute derivative of u at tn+1 given its value
    virtual shp<aVector> computeDerivative(const shp<aVector>& solnp1,
                                          double& dt) = 0;

    virtual void shiftSolutions(const shp<aVector>& sol) = 0;

    virtual shp<aVector> computeExtrapolatedSolution() = 0;

    virtual shp<aVector> combineOldSolutions() = 0;

    virtual std::vector<double> getCoefficients() const = 0;

    void dumpSolverStatistics(std::vector<SolverStatistics> statistics,
                              const double& t) const;

    void setComm(EPETRACOMM comm) {M_comm = comm; M_systemSolver.setComm(comm);}

protected:
    void initializeStatisticsFile();

    DataContainer                                       M_data;
    SystemSolver                                        M_systemSolver;
    shp<FunProvider>                                    M_funProvider;
    std::string                                         M_statisticsFile;
    EPETRACOMM                                          M_comm;
};

}

#endif // aTIMEMARCHINGALGORITHM_HPP
