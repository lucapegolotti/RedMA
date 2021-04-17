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

#ifndef ZERNIKEBASISFUNCTION_HPP
#define ZERNIKEBASISFUNCTION_HPP

#include <functional>

#include <redma/coupling_basis_functions/BasisFunctionFunctor.hpp>

namespace RedMA
{

/*! \brief A Zernike polynomial basis function.
 *
 * For the definition of the function, we refer, e.g., to
 * https://en.wikipedia.org/wiki/Zernike_polynomials.
 */
class ZernikeBasisFunction : public BasisFunctionFunctor
{
public:
    /*! \brief Constructor.
     *
     * Parameter nMax is the largest n in de set (the largest of the two indices
     * of Zernike polynomials).
     *
     * \param face the GeometricFace.
     * \param nMax Number of frequencies in the theta direction.
     */
    ZernikeBasisFunction(const GeometricFace& face,
                         unsigned int nMax);

    /*! \brief Evaluation operator.
     *
     * \param pos Position where the function has to be evaluated.
     * \return Value of the basis function.
     */
    return_Type operator()(const Vector3D& pos) override;

    /* \brief Set the index of the current basis function.
     *
     * \param index The index.
     */
    virtual void setIndex(const unsigned int& index) override;

private:
    void fillFactorials(unsigned int nMax);

    void computeOrthonormalizationCoefficient();

    unsigned int                              M_nMax;
    std::vector<int>                          M_ms;
    std::vector<int>                          M_ns;
    std::vector<std::vector<int> >            M_polyCoefs;
    std::vector<unsigned int>                 M_factorials;
    std::function<double(double)>             M_curFunction;
    double                                    M_orthoCoefficient;
    double                                    M_R;
};

}  // namespace RedMA

#endif  // FOURIERBASISFUNCTION_HPP
