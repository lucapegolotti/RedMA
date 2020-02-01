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

#ifndef LINEARSYSTEMSOLVER_HPP
#define LINEARSYSTEMSOLVER_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/array/BlockVector.hpp>
#include <redma/solver/array/BlockMatrix.hpp>
#include <redma/solver/system_solver/LinearOperatorEp.hpp>
#include <redma/solver/system_solver/InverseOperatorEp.hpp>
#include <redma/solver/system_solver/SaddlePointPreconditionerEp.hpp>

#include <redma/solver/problem/DataContainer.hpp>

#include <memory>

namespace RedMA
{

struct SolverStatistics
{
    double M_precSetupTime;
    double M_numIterations;
    double M_solveTime;
};


template<class InVectorType, class InMatrixType>
class LinearSystemSolver
{
    typedef BlockVector<InVectorType>               BV;
    typedef BlockMatrix<InMatrixType>               BM;

public:
    LinearSystemSolver(const DataContainer& datafile);

    // I don't provide  a generic implementation of this method but only
    // (template) specializations in the cpp
    void solve(const BM& matrix, const BV& rhs, BV& sol);

    void buildPreconditioner(const BM& matrix);

    SolverStatistics getSolverStatistics() const {return M_statistics;}

private:

    DataContainer                                   M_data;
    SHP(InverseOperatorEp)                          M_invOper;
    SHP(LinearOperatorEp)                           M_oper;
    SHP(PreconditionerOperatorEp)                   M_prec;
    SolverStatistics                                M_statistics;
};

}

#include "LinearSystemSolver_imp.hpp"

#endif // LINEARSYSTEMSOLVER_HPP
