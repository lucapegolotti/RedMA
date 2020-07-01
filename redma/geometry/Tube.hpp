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

#include <redma/geometry/BuildingBlock.hpp>
#include <redma/geometry/NonAffineDeformer.hpp>

namespace RedMA
{

class Tube : public BuildingBlock
{
public:
    Tube(commPtr_Type comm, std::string refinement = "coarse",
         bool verbose = false, int diameter = 1, int length = 1);

    virtual inline unsigned int expectedNumberOfChildren() override
    {
        return 1;
    }

    virtual void applyNonAffineTransformation(bool transformMesh) override;

    std::string getStringMesh(std::string refinement);

    inline double getInletRadius() const {return M_inletRadiusRef;};

    Vector3D getInletNormal() {return M_inletNormalRef;};

    double getDefLength() {return (M_outletCenterRef-M_inletCenterRef).norm();};

    std::string getOptionalParameter(unsigned int index) override;

    void resetInletOutlets() override;

    inline unsigned int getDiameter() {return M_diameter;};

    inline unsigned int getLength() {return M_length;};

    static void bendFunctionAnalytic(double& x, double& y, double& z,
                                     const double& bendAngle,
                                     const double& L);

    virtual Matrix3D computeJacobianNonAffineTransformation(const double& x,
                                                            const double& y,
                                                            const double& z) override;

private:
    void nonAffineScaling(const double& lengthRatio,
                          const double& outRadiusRatio,
                          std::shared_ptr<Transformer> transformer);

    static void scalingFunction(double& x, double& y, double& z,
                                const double& lenghtRatio,
                                const double& outRadiusRatio,
                                const double& L);

    Matrix3D computeJacobianNonAffineScaling(const double& x, const double& y, const double& z,
                                             const double& lengthRatio,
                                             const double& outRadiusRatio,
                                             const double& L);

    static double bendFunction(const double& t, const double& x, const double& y,
                               const double& z, const LifeV::ID& i,
                               const Vector3D& rotationCenter,
                               const Matrix3D& rotationMatrix);


    Matrix3D computeJacobianBend(const double& x, const double& y, const double& z,
                                 const double& bendAngle, const double& L);

    void bend(const double& bendAngle, std::shared_ptr<Transformer> transformer,
              bool transformMesh = true);

    Vector3D M_inletCenterRef;
    Vector3D M_inletNormalRef;
    Vector3D M_outletCenterRef;
    Vector3D M_outletNormalRef;

    double M_inletRadiusRef;
    double M_outletRadiusRef;

    unsigned int M_diameter;
    unsigned int M_length;

    Matrix3D M_jacobianScaling;
    Matrix3D M_jacobianBending;
};

}  // namespace RedMA

#endif  // TUBE_HPP
