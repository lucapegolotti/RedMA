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
#include <redma/solver/array/BlockVector.hpp>
#include <redma/solver/array/BlockMatrix.hpp>
#include <redma/solver/system_solver/LinearSystemSolver.hpp>

#include <sstream>
#include <memory>

namespace RedMA
{

template<class InVectorType, class InMatrixType>
class SystemSolver
{
    typedef BlockVector<InVectorType>               BV;
    typedef BlockMatrix<InMatrixType>               BM;

public:
    SystemSolver(const DataContainer& data);

    BV solve(FunctionFunctor<BV,BV> fun, FunctionFunctor<BV,BM> jac,
             BV initialGuess, int& status);

    inline std::vector<SolverStatistics> getSolverStatistics() const {return M_solverStatistics;}

    void isLinearProblem() {M_isLinearProblem = true;}

private:
    DataContainer                                       M_data;
    LinearSystemSolver<InVectorType, InMatrixType>      M_linearSystemSolver;
    std::vector<SolverStatistics>                       M_solverStatistics;
    bool                                                M_isLinearProblem;
};

}

#include "SystemSolver_imp.hpp"

#endif // SYSTEMSOLVER_HPP
