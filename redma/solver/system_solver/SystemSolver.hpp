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

#ifndef SYSTEMSOLVER_HPP
#define SYSTEMSOLVER_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/system_solver/FunctionFunctor.hpp>
#include <redma/array/BlockVector.hpp>
#include <redma/array/BlockMatrix.hpp>
#include <redma/solver/system_solver/LinearSystemSolver.hpp>

#include <sstream>
#include <memory>

namespace RedMA
{

/// Solver for a general nonlinear problem.
class SystemSolver
{
    typedef shp<aVector>               BV;
    typedef shp<aMatrix>               BM;

public:
    /*! \brief Constructor.
     *
     * \param data The dataContainer of the problem.
     */
    SystemSolver(const DataContainer& data);

    /*! \brief solve method.
     *
     * We want to solve a nonlinear equation in the form f(u) = 0.
     *
     * \param fun The FunctionFunctor of the nonlinear function.
     * \param jac The FunctionFunctor of the jacobian.
     * \param initialGuess Shared pointer to the initial guess.
     * \param status Return code; 0 if successful.
     */
    BV solve(FunctionFunctor<BV,BV> fun,
             FunctionFunctor<BV,BM> jac,
             BV initialGuess,
             int& status);

    /// Set pressure mass matrix.
    void setPressureMass(const BM& mass);

    /*! \brief Getter for the solver statistics.
     *
     * \return Vector of solver statistics (one for every iteration).
     */
    inline std::vector<SolverStatistics> getSolverStatistics() const {return M_solverStatistics;}

    /// Set M_isLinearProblem to true.
    inline void isLinearProblem() {M_isLinearProblem = true;}

    /// Setter for the MPI Communicator.
    inline void setComm(EPETRACOMM comm) {M_comm = comm;}

private:
    DataContainer                                       M_data;
    LinearSystemSolver                                  M_linearSystemSolver;
    std::vector<SolverStatistics>                       M_solverStatistics;
    bool                                                M_isLinearProblem;
    EPETRACOMM                                          M_comm;
    BM                                                  M_Mp;
};

}

#endif // SYSTEMSOLVER_HPP
