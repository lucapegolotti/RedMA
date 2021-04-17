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

/// An abstract time marching algorithm.
class aTimeMarchingAlgorithm
{
    typedef aFunctionProvider         FunProvider;
public:

    /*! \brief Constructor.
     *
     * \param datafile The DataContainer of the problem.
     */
    aTimeMarchingAlgorithm(const DataContainer& datafile);

    /*! \brief Constructor.
     *
     * \param datafile The DataContainer of the problem.
     * \param funProvider The function provider.
     */
    aTimeMarchingAlgorithm(const DataContainer& datafile,
                           shp<FunProvider> funProvider);

    /*! \brief Virtual setup function.
     *
     * \param zeroVector Shared pointer to a zero vector.
     */
    virtual void setup(const shp<aVector>& zeroVector) = 0;

    /*! \brief Virtual advance function.
     *
     * \param time The time.
     * \param dt The timestep size.
     * \param status Return code; 0 if successful.
     */
    virtual shp<aVector> advance(const double& time,
                                 double& dt,
                                 int& status) = 0;

    /*! \brief Compute derivative of a function.
     *
     * \param solnp1 Solution at time n+1.
     * \param dt Timestep size.
     * \return Shared pointer to the derivative.
     */
    virtual shp<aVector> computeDerivative(const shp<aVector>& solnp1,
                                          double& dt) = 0;

    /*! \brief Shift previous solutions given the new one.
     *
     * \param Shared pointer to the new solution.
     */
    virtual void shiftSolutions(const shp<aVector>& sol) = 0;

    /*! \brief Dump solver statistics to file.
     *
     * \param Vector of solver statistics.
     * \param time The time.
     */
    void dumpSolverStatistics(std::vector<SolverStatistics> statistics,
                              const double& time) const;

    virtual double getCoefficientExtrapolation() =0;
    virtual shp<aVector> getPreviousContribution()=0;
    virtual shp<aVector> advanceDisp(const double &dt, const shp<BlockVector> &sol)=0;
    virtual shp<aVector> getPreviousSolution()=0;



    /*! \brief Setter for the MPI Communicator.
     *
     * \param Shared pointer to the MPI Communicator.
     */

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
