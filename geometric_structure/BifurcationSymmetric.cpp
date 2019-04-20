#include "BifurcationSymmetric.hpp"

namespace RedMA
{

BifurcationSymmetric::
BifurcationSymmetric(commPtr_Type comm, bool verbose) :
  BuildingBlock(comm, verbose)
{
    M_name = "BifurcationSymmetric";

    M_datafileName = "bifurcation_symmetric_coarse_data";

    // values are referred to the reference configuration

    M_inletCenterRef[0] = 0;
    M_inletCenterRef[1] = 0;
    M_inletCenterRef[2] = 0;

    M_inletNormalRef[0] = 0;
    M_inletNormalRef[1] = -1;
    M_inletNormalRef[2] = 0;

    M_outlet1CenterRef[0] = 9;
    M_outlet1CenterRef[1] = 30;
    M_outlet1CenterRef[2] = 0;

    M_outlet1NormalRef[0] = 0.3714;
    M_outlet1NormalRef[1] = 0.9285;
    M_outlet1NormalRef[2] = 0.0;

    M_outlet2CenterRef[0] = -9;
    M_outlet2CenterRef[1] = 30;
    M_outlet2CenterRef[2] = 0.0;

    M_outlet2NormalRef[0] = -0.3714;
    M_outlet2NormalRef[1] = 0.9285;
    M_outlet2NormalRef[2] = 0.0;

    M_inletRadiusRef = 3.0;
    M_outlet1RadiusRef = 3.0;
    M_outlet2RadiusRef = 3.0;

    GeometricFace inlet(M_inletCenterRef, M_inletNormalRef, M_inletRadiusRef);
    GeometricFace outlet1(M_outlet1CenterRef, M_outlet1NormalRef, M_outlet1RadiusRef);
    GeometricFace outlet2(M_outlet2CenterRef, M_outlet2NormalRef, M_outlet2RadiusRef);

    M_inlet = inlet;
    M_outlets.push_back(outlet1);
    M_outlets.push_back(outlet2);

    // it is important to fill parametersMap right at this level because then
    // the keys will be used in the parser to check the values in the XML file
    // center of inlet
    M_parametersMap["bend"] = 0.0;
    M_parametersMap["r0"] = M_inletRadiusRef;
    M_parametersMap["r1"] = M_outlet1RadiusRef;
    M_parametersMap["r2"] = M_outlet1RadiusRef;
}

void
BifurcationSymmetric::
applyNonAffineTransformation()
{

}

}
