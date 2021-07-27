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

#ifndef INVERSEOPERATOR_HPP
#define INVERSEOPERATOR_HPP

#include <redma/RedMA.hpp>
#include <redma/problem/DataContainer.hpp>
#include <redma/solver/system_solver/LinearOperator.hpp>
#include <redma/solver/system_solver/PreconditionerOperator.hpp>

#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/linear_algebra/InvertibleOperator.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

namespace RedMA
{

/// Inverse operator.
class InverseOperator
{
    typedef LifeV::Operators::InvertibleOperator        InvertibleOperatorType;

public:
    /* \brief Constructor.
     *
     * \param data The DataContainer of the problem.
     */
    InverseOperator(const DataContainer& data);

    /* \brief Setter for the block maps.
     *
     * \param maps Shared pointer to the block maps.
     */
    void setBlockMaps(shp<BlockMaps> maps);

    /* \brief Setter for the operator.
     *
     * \param oper Shared pointer to the operator.
     */
    void setOperator(shp<LinearOperator> oper);

    /* \brief Setter for the preconditioner.
     *
     * \param prec Shared pointer to the preconditioner.
     */
    void setPreconditioner(shp<PreconditionerOperator> prec);

    /// Setter for the solver options.
    void setSolverOptions();

    /*! \brief Function to apply the inverse of the operator.
     *
     * \param rhs Shared pointer to right-hand side.
     * \param sol Shared pointer to solution.
     * \return Integer code; 0 if successful.
     */
    int invert(const shp<aVector>& rhs,
               shp<aVector>& sol);

private:
    DataContainer                                   M_data;
    shp<InvertibleOperatorType>                     M_invOper;
    shp<Teuchos::ParameterList>                     M_pListLinSolver;
    Teuchos::RCP<Teuchos::ParameterList>            M_solversOptions;
    shp<BlockMaps>                                  M_maps;
};

}

#endif // INVERSEOPERATOR_HPP
