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

#ifndef FOURIERRINGBASISFUNCTION_HPP
#define FOURIERRINGBASISFUNCTION_HPP

#include <redma/coupling_basis_functions/BasisFunctionFunctor.hpp>

namespace RedMA {

class FourierRingBasisFunction : public BasisFunctionFunctor
{
public:
    FourierRingBasisFunction(const GeometricFace& face,
                             int nFrequenciesTheta);

    return_Type operator()(const Vector3D& pos);

private:
    int                         M_nFrequenciesTheta;
    std::vector<double>         M_thetaFreq;
    std::vector<double>         M_thetaPhase;
    std::vector<unsigned int>   M_auxIndicesTheta;

    double                      M_sigmaRadial;
};

} // Namespace RedMA


#endif // FOURIERRINGBASISFUNCTION_HPP
