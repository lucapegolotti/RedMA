#include "Tube.hpp"

namespace RedMA
{

Tube::
Tube(commPtr_Type comm, bool verbose) :
  BuildingBlock(comm, verbose)
{
    M_name = "Tube";

    M_datafileName = "tube_coarse_data";

    // center of inlet (reference configuration)
    M_inletCenterRef[0] = 0.0;
    M_inletCenterRef[1] = 0.0;
    M_inletCenterRef[2] = 0.0;

    // center of outlet (reference configuration)
    M_outletCenterRef[0] = 0.0;
    M_outletCenterRef[1] = 0.0;
    M_outletCenterRef[2] = 15.0;

    // normal of inlet (reference configuration)
    M_inletNormalRef[0] = 0.0;
    M_inletNormalRef[1] = 0.0;
    M_inletNormalRef[2] = -1.0;

    // outlet of outlet (reference configuration)
    M_outletNormalRef[0] = 0.0;
    M_outletNormalRef[1] = 0.0;
    M_outletNormalRef[2] = 1.0;

    M_inletRadiusRef = 3.0;
    M_outletRadiusRef = 3.0;

    GeometricFace inlet(M_inletCenterRef, M_inletNormalRef, M_inletRadiusRef);
    GeometricFace outlet(M_outletCenterRef, M_outletNormalRef, M_outletRadiusRef);

    M_inlet = inlet;
    M_outlets.push_back(outlet);

    // it is important to fill parametersMap right at this level because then
    // the keys will be used in the parser to check the values in the XML file
    // center of inlet
    M_parametersMap["bend"] = 0.0;
    M_parametersMap["L_ratio"] = 1.0;
    M_parametersMap["Rout_ratio"] = 1.0;
}

void
Tube::
applyNonAffineTransformation()
{
    LifeV::MeshUtility::MeshTransformer<mesh_Type> transformer(*M_mesh);
    nonAffineScaling(M_parametersMap["L_ratio"], M_parametersMap["Rout_ratio"],
                     transformer);
}

void
Tube::
nonAffineScaling(const double lengthRatio, const double radiusRatio,
                 Tube::Transformer& transformer)
{
    auto foo = std::bind(scalingFunction, std::placeholders::_1,
                         std::placeholders::_2, std::placeholders::_3,
                         lengthRatio, radiusRatio);
    transformer.transformMesh(foo);

    M_outlets[0].M_radius = M_outlets[0].M_radius * radiusRatio;
    M_outlets[0].M_center[2] = M_outlets[0].M_center[2] * lengthRatio;
}

void
Tube::
scalingFunction(double& x, double& y, double& z,
                const double lenghtRatio, const double outRadiusRatio)
{
    // 15 is the length of the tube
    double curRatio = 1. - (1. - outRadiusRatio) * z / 15.0;
    z = z * lenghtRatio;
    x = x * curRatio;
    y = y * curRatio;
}

}
