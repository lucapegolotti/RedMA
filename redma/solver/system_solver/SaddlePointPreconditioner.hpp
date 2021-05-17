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
#include <lifev/navier_stokes_blocks/solver/aPmmOperator.hpp>
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
    typedef shp<BlockVector>                                         BV;
    typedef shp<BlockMatrix>                                         BM;
    typedef LifeV::Operators::NavierStokesPreconditionerOperator     NSPrec;
    typedef LifeV::Operators::NavierStokesOperator                   NSOp;
    typedef LifeV::Operators::InvertibleOperator                     InvOp;
    typedef LifeV::Operators::ApproximatedInvertibleRowMatrix        ApproxInv;

public:
    SaddlePointPreconditioner(const DataContainer& data,
                              const BM& matrix,
                              const BM& pressureMass = nullptr);

    virtual int ApplyInverse(const super::vector_Type& X,
                             super::vector_Type& Y) const override;

    void allocateInnerPreconditioners(const BM& primalMatrix);

    void allocateInverseSolvers(const BM& primalMatrix);

    void allocateApproximatedInverses(const BM& primalMatrix);

    void setSolverOptions();

    void computeSchurComplement(const BM& A, const BM& BT,
                                const BM& B, const BM& C);

    void setup(const BM& matrix,
               const BM& pressureMass = nullptr,
               bool doComputeSchurComplement = true);

private:
    void computeAm1BT(const BM& A, const BM& BT);

    shp<aMatrix> computeSingleAm1BT(const BM& A,
                                    const BM& BT,
                                    const unsigned int& index);

    void solveEveryPrimalBlock(const VECTOREPETRA& X, VECTOREPETRA &Y) const;

    void applyEveryB(const VECTOREPETRA& X, VECTOREPETRA &Y) const;

    void applyEveryAm1BT(const VECTOREPETRA& X, VECTOREPETRA &Y) const;

    void findSmallBlocks(const BM& primalMatrix);

    DataContainer                                        M_data;
    BM                                                   M_matrixCollapsed;
    BM                                                   M_S;
    shp<BlockMaps>                                       M_maps;
    std::vector<shp<NSPrec>>                             M_innerPreconditioners;
    std::vector<shp<InvOp>>                              M_invOperators;
    shp<Teuchos::ParameterList>                          M_pListLinSolver;
    Teuchos::RCP<Teuchos::ParameterList>                 M_solversOptionsInner;
    BM                                                   M_Am1BT;
    std::string                                          M_innerPrecType;
    std::string                                          M_approxSchurType;
    shp<ApproxInv>                                       M_approximatedSchurInverse;
    std::vector<shp<ApproxInv>>                          M_approximatedInverses;
    shp<MAPEPETRA>                                       M_primalMap;
    shp<MAPEPETRA>                                       M_dualMap;
    shp<MAPEPETRA>                                       M_monolithicMap;
    std::vector<shp<MAPEPETRA>>                          M_rangeMaps;
    std::vector<shp<MAPEPETRA>>                          M_domainMaps;
    unsigned int                                         M_nPrimalBlocks;
    unsigned int                                         M_nDualBlocks;
    std::vector<shp<Epetra_SerialDenseSolver>>           M_solversAsDense;
    int                                                  M_thresholdSizeExactSolve;
    std::vector<bool>                                    M_isSmallBlock;
};

}

#endif // SADDLEPOINTPRECONDITIONER_HPP
