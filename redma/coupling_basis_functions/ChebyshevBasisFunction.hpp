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

#ifndef CHEBYSHEVBASISFUNCTION_HPP
#define CHEBYSHEVBASISFUNCTION_HPP

#include <functional>

#include <redma/coupling_basis_functions/BasisFunctionFunctor.hpp>

namespace RedMA
{

class ChebyshevBasisFunction : public BasisFunctionFunctor
{
public:
    ChebyshevBasisFunction(const GeometricFace& face,
                           unsigned int nMax);

    return_Type operator()(const Vector3D& pos) override;

private:
    double chebyshevU(const double& x, const unsigned int& n);

    unsigned int                              M_nMax;
    std::vector<int>                          M_ks;
    std::vector<int>                          M_ns;
    double                                    M_R;
    // versor to compute the normal
    Vector3D                                  M_e;
    const double                              M_sqrtPIm1;
};

}  // namespace RedMA

#endif  // CHEBYSHEVBASISFUNCTION_HPP
