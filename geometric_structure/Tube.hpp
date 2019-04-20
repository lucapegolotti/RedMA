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

namespace RedMA
{

class Tube : public BuildingBlock
{
    typedef LifeV::MeshUtility::MeshTransformer<mesh_Type> Transformer;
public:
    Tube(commPtr_Type comm, bool verbose = false);

    virtual inline unsigned int expectedNumberOfChildren()
    {
        return 1;
    }

    virtual void applyNonAffineTransformation();

private:
    void nonAffineScaling(const double lengthRatio,
                          const double outRadiusRatio,
                          Transformer& transformer);

    static void scalingFunction(double& x, double& y, double& z,
                                const double lenghtRatio,
                                const double outRadiusRatio);

    Vector3D M_inletCenterRef;
    Vector3D M_inletNormalRef;
    Vector3D M_outletCenterRef;
    Vector3D M_outletNormalRef;

    double M_inletRadiusRef;
    double M_outletRadiusRef;
};

}  // namespace RedMA

#endif  // TUBE_HPP
