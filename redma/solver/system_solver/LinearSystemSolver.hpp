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
#include <redma/array/BlockVector.hpp>
#include <redma/array/BlockMatrix.hpp>
#include <redma/solver/system_solver/LinearOperator.hpp>
#include <redma/solver/system_solver/InverseOperator.hpp>
#include <redma/solver/system_solver/SaddlePointPreconditioner.hpp>

#include <redma/problem/DataContainer.hpp>

#include <memory>

namespace RedMA
{

struct SolverStatistics
{
    double M_precSetupTime;
    double M_numIterations;
    double M_solveTime;
};

class LinearSystemSolver
{
    typedef shp<aVector>               BV;
    typedef shp<aMatrix>               BM;

public:
    LinearSystemSolver(const DataContainer& datafile);

    void solve(const BM& matrix, const BV& rhs, BV& sol);

    void buildPreconditioner(const BM& matrix);

    SolverStatistics getSolverStatistics() const {return M_statistics;}

    void setComm(EPETRACOMM comm) {M_comm = comm;}

private:

    // only required for dense computation
    void computeSchurComplementDense(const BM& matrix);

    void solveDense(const BM& matrix, const BV& rhs, BV& sol);

    void convertVectorType(const shp<BlockMatrix>& matrix,
                           const shp<DenseVector>& vector,
                           shp<BlockVector>& targetVector);

    // only required for dense computation
    std::vector<shp<Epetra_SerialDenseSolver>>      M_solversAs;
    Epetra_SerialDenseSolver                        M_schurSolver;
    std::vector<shp<DenseMatrix>>                   M_collapsedAs;
    shp<DenseMatrix>                                M_schurComplementColl;

    DataContainer                                   M_data;
    shp<InverseOperator>                            M_invOper;
    shp<LinearOperator>                             M_oper;
    shp<PreconditionerOperator>                     M_prec;
    shp<BlockMaps>                                  M_maps;

    SolverStatistics                                M_statistics;
    unsigned                                        M_numSolves;

    EPETRACOMM                                      M_comm;
};

}

#endif // LINEARSYSTEMSOLVER_HPP
