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

#ifndef TUBE_HPP
#define TUBE_HPP

#include <BuildingBlock.hpp>
#include <NonAffineDeformer.hpp>

namespace RedMA
{

class Tube : public BuildingBlock
{
public:
    Tube(commPtr_Type comm, std::string refinement = "coarse",
         bool verbose = false, int diameter = 1, int length = 1);

    virtual inline unsigned int expectedNumberOfChildren()
    {
        return 1;
    }

    virtual void applyNonAffineTransformation();

    std::string getStringMesh(std::string refinement);

    inline double getInletRadius() const {return M_inletRadiusRef;};

    Vector3D getInletNormal() {return M_inletNormalRef;};

private:
    void nonAffineScaling(const double& lengthRatio,
                          const double& outRadiusRatio,
                          Transformer& transformer);

    static void scalingFunction(double& x, double& y, double& z,
                                const double lenghtRatio,
                                const double outRadiusRatio);

    static double bendFunction(const double& t, const double& x, const double& y,
                               const double& z, const LifeV::ID& i,
                               const Vector3D& rotationCenter,
                               const Matrix3D& rotationMatrix);

    static void bendFunctionAnalytic(double& x, double& y, double& z,
                                     const double& bendAngle,
                                     const GeometricFace& outlet);

    void bend(const double& bendAngle, Transformer& transformer);

    Vector3D M_inletCenterRef;
    Vector3D M_inletNormalRef;
    Vector3D M_outletCenterRef;
    Vector3D M_outletNormalRef;

    double M_inletRadiusRef;
    double M_outletRadiusRef;
};

}  // namespace RedMA

#endif  // TUBE_HPP
