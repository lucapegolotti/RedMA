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

#ifndef INVERSEOPERATOREP_HPP
#define INVERSEOPERATOREP_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/problem/DataContainer.hpp>
#include <redma/solver/system_solver/LinearOperatorEp.hpp>
#include <redma/solver/system_solver/PreconditionerOperatorEp.hpp>

#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/linear_algebra/InvertibleOperator.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

namespace RedMA
{

class InverseOperatorEp
{
    typedef LifeV::Operators::InvertibleOperator        InvertibleOperatorType;

public:
    InverseOperatorEp(const DataContainer& data);

    void setOperator(SHP(LinearOperatorEp) oper);

    void setPreconditioner(SHP(PreconditionerOperatorEp) prec);

    void setSolverOptions();

    void invert(const BlockVector<BlockVector<VectorEp>>& rhs,
                BlockVector<BlockVector<VectorEp>>& sol);

private:
    DataContainer                                   M_data;
    SHP(InvertibleOperatorType)                     M_invOper;
    SHP(Teuchos::ParameterList)                     M_pListLinSolver;
    Teuchos::RCP<Teuchos::ParameterList>            M_solversOptions;
    SHP(BlockMaps<BlockMatrix<MatrixEp>>)           M_maps;
};

}

#endif // INVERSEOPERATOREP_HPP
