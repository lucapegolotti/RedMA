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

#include <redma/geometry/BuildingBlock.hpp>
#include <redma/geometry/NonAffineDeformer.hpp>

namespace RedMA
{

class BifurcationSymmetric : public BuildingBlock
{
public:
    BifurcationSymmetric(commPtr_Type comm, std::string refinement = "coarse",
                         bool verbose = false, int angle = 50, bool randomizable = true);

    virtual inline unsigned int expectedNumberOfChildren() override
    {
        return 2;
    }

    virtual void applyNonAffineTransformation(bool transformMesh = true) override;

    inline Vector3D getCenter() const {return M_center;};

    inline Vector3D getInletNormal() const {return M_inletNormalRef;};

    inline Vector3D getTransverse() const {return M_transverse;};

    inline double getInletRadius() const {return M_inletRadiusRef;};

    void resetInletOutlets() override;

    virtual std::string getOptionalParameter(unsigned int index) override;

    virtual Matrix3D computeJacobianNonAffineTransformation(const double& x,
                                                            const double& y,
                                                            const double& z) override {return M_identity3D;};

private:
    static double outletMapFunction(const double& t, const double& x,
                               const double& y, const double& z,
                               const LifeV::ID& i,
                               const GeometricFace& targetFace,
                               const Vector3D& desiredCenter,
                               const Matrix3D& rotationMatrix);

    void bend(const double& out1_alphax,
              const double& out1_alphay,
              const double& out1_alphaz,
              const double& out2_alphax,
              const double& out2_alphay,
              const double& out2_alphaz,
              std::shared_ptr<Transformer> transformer,
              bool transformMesh = true);

    void rotateGeometricFace(const GeometricFace& face, Vector3D& rotatedCenter,
                             Vector3D& rotatedNormal,
                             const Matrix3D& rotationMatrix,
                             const Vector3D& rotationCenter);

    void computeCenter();

    std::string getStringMesh(std::string refinement);

    Vector3D M_inletCenterRef;
    Vector3D M_inletNormalRef;
    Vector3D M_outlet1CenterRef;
    Vector3D M_outlet1NormalRef;
    Vector3D M_outlet2CenterRef;
    Vector3D M_outlet2NormalRef;

    Vector3D M_center;

    Vector3D M_transverse;

    Matrix3D M_identity3D;

    double M_inletRadiusRef;
    double M_outlet1RadiusRef;
    double M_outlet2RadiusRef;
    int    M_angle;
};

}  // namespace RedMA

#endif  // BIFURCATIONSYMMETRIC_HPP
