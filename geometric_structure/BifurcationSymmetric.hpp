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

#ifndef BIFURCATIONSYMMETRIC_HPP
#define BIFURCATIONSYMMETRIC_HPP

#include "BuildingBlock.hpp"

namespace RedMA
{

class BifurcationSymmetric : public BuildingBlock
{
public:
    BifurcationSymmetric(commPtr_Type comm, bool verbose = false);

    virtual inline unsigned int expectedNumberOfChildren()
    {
        return 2;
    }

    virtual void applyNonAffineTransformation();

private:
    static double outletMapFunction(const double& t, const double& x,
                               const double& y, const double& z,
                               const LifeV::ID& i,
                               const GeometricFace& targetFace,
                               const Vector3D& desiredCenter,
                               const Matrix3D& rotationMatrix);

    void bend(const double& out1_alpha_plane, const double& out1_alphaz,
              const double& out2_alpha_plane, const double& out2_alphaz,
              Transformer& transformer);

    void rotateGeometricFace(const GeometricFace& face, Vector3D& rotatedCenter,
                             Vector3D& rotatedNormal,
                             const Matrix3D& rotationMatrix,
                             const Vector3D& rotationCenter);

    Vector3D M_inletCenterRef;
    Vector3D M_inletNormalRef;
    Vector3D M_outlet1CenterRef;
    Vector3D M_outlet1NormalRef;
    Vector3D M_outlet2CenterRef;
    Vector3D M_outlet2NormalRef;

    double M_inletRadiusRef;
    double M_outlet1RadiusRef;
    double M_outlet2RadiusRef;
};

}  // namespace RedMA

#endif  // BIFURCATIONSYMMETRIC_HPP
