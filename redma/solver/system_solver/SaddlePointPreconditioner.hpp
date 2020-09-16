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

#ifndef SADDLEPOINTPRECONDITIONER_HPP
#define SADDLEPOINTPRECONDITIONER_HPP

#include <redma/solver/system_solver/PreconditionerOperator.hpp>
#include <redma/problem/DataContainer.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/utils/PrintLog.hpp>

#include <lifev/navier_stokes_blocks/solver/NavierStokesOperator.hpp>
#include <lifev/navier_stokes_blocks/solver/aSIMPLEOperator.hpp>
#include <lifev/core/linear_algebra/LinearOperatorAlgebra.hpp>
#include <lifev/core/linear_algebra/BlockEpetra_MultiVector.hpp>
#include <lifev/core/linear_algebra/BlockEpetra_Map.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

namespace RedMA
{

class SaddlePointPreconditioner : public PreconditionerOperator
{
    typedef LifeV::Operators::LinearOperatorAlgebra                  super;
    typedef SHP(aVector)                                             BV;
    typedef SHP(aMatrix)                                             BM;
    typedef LifeV::Operators::NavierStokesPreconditionerOperator     NSPrec;
    typedef LifeV::Operators::NavierStokesOperator                   NSOp;
    typedef LifeV::Operators::InvertibleOperator                     InvOp;
    typedef LifeV::Operators::ApproximatedInvertibleRowMatrix        ApproxInv;

public:
    SaddlePointPreconditioner(const DataContainer& data, const BM& matrix);

    virtual int ApplyInverse(const super::vector_Type& X,
                             super::vector_Type& Y) const override;

    void allocateInnerPreconditioners(const BM& primalMatrix);

    void allocateInverseSolvers(const BM& primalMatrix);

    void setSolverOptions();


    void computeSchurComplement(const BM& A, const BM& BT,
                                const BM& B, const BM& C);

private:
    void computeAm1BT(const BM& A, const BM& BT);

    SHP(aMatrix) computeSingleAm1BT(const BM& A,
                                    const BM& BT,
                                    const unsigned int& index);

    void solveEveryPrimalBlock(const VECTOREPETRA& X, VECTOREPETRA &Y) const;

    void applyEveryB(const VECTOREPETRA& X, VECTOREPETRA &Y) const;

    void applyEveryAm1BT(const VECTOREPETRA& X, VECTOREPETRA &Y) const;

    DataContainer                                        M_data;
    BM                                                   M_matrix;
    BM                                                   M_matrixCollapsed;
    BM                                                   M_S;
    std::vector<SHP(NSPrec)>                             M_innerPreconditioners;
    std::vector<SHP(InvOp)>                              M_invOperators;
    SHP(Teuchos::ParameterList)                          M_pListLinSolver;
    Teuchos::RCP<Teuchos::ParameterList>                 M_solversOptionsInner;
    BM                                                   M_Am1BT;
    std::string                                          M_innerPrecType;
    std::string                                          M_approxSchurType;
    SHP(ApproxInv)                                       M_approximatedSchurInverse;
    SHP(MAPEPETRA)                                       M_primalMap;
    SHP(MAPEPETRA)                                       M_dualMap;
    SHP(MAPEPETRA)                                       M_monolithicMap;
    std::vector<SHP(MAPEPETRA)>                          M_rangeMaps;
    std::vector<SHP(MAPEPETRA)>                          M_domainMaps;
    unsigned int                                         M_nPrimalBlocks;
    unsigned int                                         M_nDualBlocks;
};

}

#endif // SADDLEPOINTPRECONDITIONER_HPP
