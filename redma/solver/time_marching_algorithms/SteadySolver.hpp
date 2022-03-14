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

#ifndef STEADYSOLVER_HPP
#define STEADYSOLVER_HPP

#include <redma/RedMA.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/array/BlockVector.hpp>
#include <redma/solver/time_marching_algorithms/aFunctionProvider.hpp>
#include <redma/solver/system_solver/FunctionFunctor.hpp>
#include <redma/solver/system_solver/SystemSolver.hpp>

namespace RedMA {

class SteadySolver {

    typedef aFunctionProvider      FunProvider;
    typedef shp<aVector>           BV;
    typedef shp<aMatrix>           BM;

public:
    /*! \brief Constructor.
     *
     * \param datafile The DataContainer of the problem.
     */
    SteadySolver(const DataContainer& data);

    /*! \brief Constructor.
     *
     * \param datafile The DataContainer of the problem.
     * \param funProvider The function provider.
     */
    SteadySolver(const DataContainer& data,
                 shp<FunProvider> funProvider);

    /*! \brief Set the initial guess from file
     *
     * \param ICpath path to the folder storing the IC. If empty, the IC is set to 0.
     */
    void setInitialGuess(const std::string& ICpath = "");

    /*! \brief Set the initial guess from a vector
     *
     * \param IC initial guess
     */
    void setInitialGuess(const BV IC);

    /*! \brief Get the initial guess
     *
     */
    BV getInitialGuess() const {return M_initialGuess;};

    /*! \brief Solve function.
     *
     * \param status Return code; 0 if successful.
     */
    shp<aVector> solve(int& status);

    /*! \brief Dump solver statistics to file.
     *
     * \param Vector of solver statistics.
     * \param time The time.
     */
    void dumpSolverStatistics(std::vector<SolverStatistics> statistics) const;

    /*! \brief Setter for the MPI Communicator.
     *
     * \param Shared pointer to the MPI Communicator.
     */
    void setComm(EPETRACOMM comm) {M_comm = comm; M_systemSolver.setComm(comm);}

    /// Set the solver as linear
    inline void setLinearSolver() {this->M_systemSolver.isLinearProblem();}

protected:
    void initializeStatisticsFile();

    DataContainer                                       M_data;
    SystemSolver                                        M_systemSolver;
    shp<FunProvider>                                    M_funProvider;
    BV                                                  M_initialGuess;
    std::string                                         M_statisticsFile;
    EPETRACOMM                                          M_comm;

};

}


#endif
