#include "Bypass.hpp"

namespace RedMA
{

Bypass::
Bypass(EPETRACOMM comm, std::string name, bool verbose) :
  BuildingBlock(comm, "coarse", verbose)
{
    M_name = name;
    M_datafileName = "data_mesh";
    M_meshName = "others/bypass_coarse_fluid.mesh";

    // center of inlet (reference configuration)
    M_inletCenterRef1[0] = -7.06006787;
    M_inletCenterRef1[1] = 1.63241487;
    M_inletCenterRef1[2] = 45.36727996;

    // center of outlet (reference configuration)
    M_inletCenterRef2[0] = -8.51731429;
    M_inletCenterRef2[1] = 2.71285607;
    M_inletCenterRef2[2] = 45.23909379;

    // center of outlet (reference configuration)
    M_outletCenterRef[0] = -8.49107192;
    M_outletCenterRef[1] = 2.23961129;
    M_outletCenterRef[2] = 40.33716354;

    // normal of inlet (reference configuration)
    M_inletNormalRef1[0] = 0.12364816;
    M_inletNormalRef1[1] = -0.21329954;
    M_inletNormalRef1[2] = 0.96913076;

    // outlet of outlet (reference configuration)
    M_inletNormalRef2[0] = -0.13280704;
    M_inletNormalRef2[1] = 0.30002537;
    M_inletNormalRef2[2] = 0.94464124;

    // outlet of outlet (reference configuration)
    M_outletNormalRef[0] = 0.27216076;
    M_outletNormalRef[1] = 0.1600235;
    M_outletNormalRef[2] = 0.94885246;

    M_inletRadiusRef1 = 0.33609344929489127;
    M_inletRadiusRef2 = 0.2661432021360946;
    M_outletRadiusRef = 0.24480087734232522;


    M_wallFlag = 200;

    resetInletOutlets();
}

void
Bypass::
resetInletOutlets()
{
    GeometricFace inlet1(M_inletCenterRef1, M_inletNormalRef1, M_inletRadiusRef1, 2, 20);
    GeometricFace inlet2(M_inletCenterRef2, M_inletNormalRef2, M_inletRadiusRef2, 3, 30);
    GeometricFace outlet(M_outletCenterRef, M_outletNormalRef, M_outletRadiusRef, 4, 40);

    M_inlets.clear();
    M_inlets.push_back(inlet1);
    M_inlets.push_back(inlet2);
    M_outlets.clear();
    M_outlets.push_back(outlet);
}


std::string
Bypass::
getOptionalParameter(unsigned int index)
{
}

void
Bypass::
applyNonAffineTransformation(bool transformMesh)
{
}

}
