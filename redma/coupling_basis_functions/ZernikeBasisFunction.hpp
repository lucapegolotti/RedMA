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

class ZernikeBasisFunction : public BasisFunctionFunctor
{
public:
    ZernikeBasisFunction(const GeometricFace& face,
                         unsigned int nMax);

    return_Type operator()(const Vector3D& pos) override;

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
