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

#ifndef SADDLEPOINTPRECONDITIONEREP_HPP
#define SADDLEPOINTPRECONDITIONEREP_HPP

#include <redma/solver/system_solver/PreconditionerOperatorEp.hpp>
#include <redma/solver/problem/DataContainer.hpp>
#include <redma/utils/Exception.hpp>

#include <lifev/navier_stokes_blocks/solver/NavierStokesOperator.hpp>
#include <lifev/navier_stokes_blocks/solver/aSIMPLEOperator.hpp>
#include <lifev/core/linear_algebra/LinearOperatorAlgebra.hpp>
#include <lifev/core/linear_algebra/BlockEpetra_MultiVector.hpp>
#include <lifev/core/linear_algebra/BlockEpetra_Map.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

namespace RedMA
{

class SaddlePointPreconditionerEp : public PreconditionerOperatorEp
{
    typedef LifeV::Operators::LinearOperatorAlgebra                  super;
    typedef BlockVector<BlockVector<VectorEp>>                       BV;
    typedef BlockMatrix<BlockMatrix<MatrixEp>>                       BM;
    typedef LifeV::Operators::NavierStokesPreconditionerOperator     NSPrec;

public:
    SaddlePointPreconditionerEp(const DataContainer& data, const BM& matrix);

    virtual int ApplyInverse(const super::vector_Type& X,
                             super::vector_Type& Y) const override;

    void allocateInnerPreconditioners(const BM& primalMatrix, std::string precType);

    void setSolverOptions();

private:
    DataContainer                                       M_data;
    BM                                                  M_matrix;
    std::vector<SHP(NSPrec)>                            M_innerPreconditioners;
    SHP(Teuchos::ParameterList)                         M_pListLinSolver;
    Teuchos::RCP<Teuchos::ParameterList>                M_solversOptionsInner;
};

}

#endif // SADDLEPOINTPRECONDITIONEREP_HPP
