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

#ifndef AORTABIFURCATION0_HPP
#define AORTABIFURCATION0_HPP

#include <redma/geometry/BuildingBlock.hpp>
#include <redma/geometry/NonAffineDeformer.hpp>

namespace RedMA
{

class AortaBifurcation0 : public BuildingBlock
{
public:
    AortaBifurcation0(commPtr_Type comm, std::string refinement,
                      std::string name = "aorta", bool verbose = false);

    virtual inline unsigned int expectedNumberOfChildren() override
    {
        return 4;
    }

    virtual void applyNonAffineTransformation(bool transformMesh) override;

    std::string getOptionalParameter(unsigned int index) override;

    void resetInletOutlets() override;

    virtual Matrix3D computeJacobianNonAffineTransformation(const double& x,
                                                            const double& y,
                                                            const double& z) override {};

private:

    Vector3D M_inletCenterRef;
    Vector3D M_inletNormalRef;
    Vector3D M_outletCenterRef1;
    Vector3D M_outletNormalRef1;
    Vector3D M_outletCenterRef2;
    Vector3D M_outletNormalRef2;
    Vector3D M_outletCenterRef3;
    Vector3D M_outletNormalRef3;
    Vector3D M_outletCenterRef4;
    Vector3D M_outletNormalRef4;

    double M_inletRadiusRef;
    double M_outletRadiusRef1;
    double M_outletRadiusRef2;
    double M_outletRadiusRef3;
    double M_outletRadiusRef4;
};

}  // namespace RedMA

#endif  // AORTA_HPP
