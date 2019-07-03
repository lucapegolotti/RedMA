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

#ifndef BASISFUNCTIONFUNCTOR_HPP
#define BASISFUNCTIONFUNCTOR_HPP

#include <functional>

#include <lifev/core/array/VectorSmall.hpp>

#include <BuildingBlock.hpp>

namespace RedMA
{

class BasisFunctionFunctor
{
public:
    typedef double                                            return_Type;
    typedef std::function<double(const double, const double)> Function;
    typedef LifeV::VectorSmall<3>                             Vector3D;

    BasisFunctionFunctor(const GeometricFace& face);

    virtual return_Type operator()(const Vector3D& pos) = 0;

    virtual void setIndex(const unsigned int& index);

    unsigned int getNumBasisFunctions() const;

protected:
    GeometricFace           M_face;
    unsigned int            M_index;
    unsigned int            M_nBasisFunctions;

};

}  // namespace RedMA

#endif  // BASISFUNCTIONFUNCTOR_HPP
