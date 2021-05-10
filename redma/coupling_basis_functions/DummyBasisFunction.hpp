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

// these are taken from https://arxiv.org/pdf/1701.02709.pdf (formula 3.8)

#ifndef DUMMYBASISFUNCTION_HPP
#define DUMMYBASISFUNCTION_HPP

#include <functional>

#include <redma/coupling_basis_functions/BasisFunctionFunctor.hpp>
#include <redma/utils/Exception.hpp>

namespace RedMA
{

/*! \brief A dummy basis function used for testing.
 *
 * The function returns 0.
 */
class DummyBasisFunction : public BasisFunctionFunctor
{
public:
    /*! \brief Constructor taking a GeometricFace.
     *
     * \param face The GeometricFace where the basis functions are defined.
     * \param type The name of the basis functions.
     */
    DummyBasisFunction(const GeometricFace& face,
                       std::string type);

    /*! \brief Evaluation operator.
     *
     * \param pos Position where the function has to be evaluated.
     * \return Value of the basis function, i.e., 0.
     */
    return_Type operator()(const Vector3D& pos) override;
};

}  // namespace RedMA

#endif  // DUMMYBASISFUNCTION_HPP
