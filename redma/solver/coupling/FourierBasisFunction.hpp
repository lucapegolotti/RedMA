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

#include <redma/solver/coupling/BasisFunctionFunctor.hpp>

namespace RedMA
{

class FourierBasisFunction : public BasisFunctionFunctor
{
public:
    FourierBasisFunction(const GeometricFace& face,
                         unsigned int nFrequenciesTheta,
                         unsigned int nFrequenciesRadial);

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
    // versor to compute the normal
    Vector3D                    M_e;
};

}  // namespace RedMA

#endif  // FOURIERBASISFUNCTION_HPP
