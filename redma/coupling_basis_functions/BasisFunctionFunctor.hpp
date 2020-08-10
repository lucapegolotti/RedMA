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

#include <redma/geometry/BuildingBlock.hpp>

namespace RedMA
{

class BasisFunctionFunctor
{
public:
    typedef double                                            return_Type;
    typedef LifeV::VectorSmall<3>                             Vector3D;
    typedef std::function<double(double const&,
                                 double const&,
                                 double const&,
                                 double const&,
                                 unsigned int const& )>       Function;

    BasisFunctionFunctor(const GeometricFace& face);

    virtual return_Type operator()(const Vector3D& pos) = 0;

    Function function();

    return_Type evaluateOperator(const double& t, const double& x, const double& y,
                                 const double& z, unsigned int const& index);

    virtual void setIndex(const unsigned int& index);

    unsigned int getNumBasisFunctions() const;

    std::string getType() const {return M_type;};

protected:
    void getLocalXAndY(const Vector3D& pos, double& x, double& y);

    void getThetaAndRadius(const Vector3D& pos, double& theta, double& radius);

    GeometricFace           M_face;
    unsigned int            M_index;
    unsigned int            M_nBasisFunctions;
    // versor to compute the normal
    Vector3D                M_e;
    Vector3D                M_eOrth;
    std::string             M_type;
};

}  // namespace RedMA

#endif  // BASISFUNCTIONFUNCTOR_HPP
