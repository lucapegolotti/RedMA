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

#ifndef FOURIERBASISFUNCTION_HPP
#define FOURIERBASISFUNCTION_HPP

#include <redma/coupling_basis_functions/BasisFunctionFunctor.hpp>

namespace RedMA
{

/*! \brief A basis functions obtain as a tensor product of Fourier bfs in the
 * \f$\theta\f$ and radial direction \f$r\f$.
 *
 * In the expressions that follow, \f$R\f$ is the radius of the face,
 * \f$N_{\theta}\f$ is nFrequenciesTheta (set in the constructor), and
 * \f$N_{r}\f$ is nFrequenciesRadial (set in the constructor).
 * The basis functions are written
 * \f{eqnarray*}{
 *        \xi_j(\theta,r) = \sin(\omega_{\theta,j} \theta + \alpha_{\theta,j}) \cos(\omega_{r,j} r + \alpha_{r,j}), \\
 * \f}
 * where \f$\omega_{\theta,0} = 0\f$,
 * \f$\omega_{\theta,2k-1} = 2k\f$,
 * \f$\omega_{\theta,2k} = 2k\f$,
 * \f$\omega_{\theta,0} = pi/2\f$,
 * \f$\alpha_{\theta,2k-1} = 0\f$,
 * \f$\alpha_{\theta,2k} = \pi/2\f$,
 * for \f$k = 1,\ldots,N_{\theta}\f$,
 * \f$\omega_{r,k} = 2 k R\f$, \f$\alpha_{r,k} = 0\f$, for \f$k = 0,\ldots,N_{r}\f$.
 *
 */
class FourierBasisFunction : public BasisFunctionFunctor
{
public:
    /*! \brief Constructor taking a GeometricFace.
     *
     * \param face the GeometricFace.
     * \param nFrequenciesTheta Number of frequencies in the theta direction.
     * \param nFrequenciesRadial Number of frequencies in the radial direction.
     */
    FourierBasisFunction(const GeometricFace& face,
                         unsigned int nFrequenciesTheta,
                         unsigned int nFrequenciesRadial);

    /*! \brief Evaluation operator.
     *
     * \param pos Position where the function has to be evaluated.
     * \return Value of the basis function.
     */
    return_Type operator()(const Vector3D& pos);

private:
    unsigned int                M_nFrequenciesTheta;
    unsigned int                M_nFrequenciesRadial;
    std::vector<double>         M_thetaFreq;
    std::vector<double>         M_radialFreq;
    std::vector<double>         M_thetaPhase;
    std::vector<double>         M_radialPhase;
    std::vector<unsigned int>   M_auxIndicesTheta;
    std::vector<unsigned int>   M_auxIndicesRadial;
};

}  // namespace RedMA

#endif  // FOURIERBASISFUNCTION_HPP
