#include "Aorta.hpp"

namespace RedMA
{

Aorta::
Aorta(commPtr_Type comm, std::string name, bool verbose) :
  BuildingBlock(comm, "normal", verbose)
{
    M_name = name;
    M_datafileName = "data_mesh";
    M_meshName = "others/aorta.mesh";

    // center of inlet (reference configuration)
    M_inletCenterRef[0] = -1.966975;
    M_inletCenterRef[1] = -1.624933;
    M_inletCenterRef[2] = 12.923909;

    // center of outlet (reference configuration)
    M_outletCenterRef1[0] = -5.716039;
    M_outletCenterRef1[1] = 5.621785;
    M_outletCenterRef1[2] = -20.529217;

    // center of outlet (reference configuration)
    M_outletCenterRef2[0] = 8.266322;
    M_outletCenterRef2[1] = 4.993420;
    M_outletCenterRef2[2] = -19.989317;

    // normal of inlet (reference configuration)
    M_inletNormalRef[0] = -0.449215;
    M_inletNormalRef[1] = -0.382102;
    M_inletNormalRef[2] = 0.807591;

    // outlet of outlet (reference configuration)
    M_outletNormalRef1[0] = -0.318503;
    M_outletNormalRef1[1] = 0.325307;
    M_outletNormalRef1[2] = -0.890355;

    // outlet of outlet (reference configuration)
    M_outletNormalRef2[0] = 0.272467;
    M_outletNormalRef2[1] = 0.387918;
    M_outletNormalRef2[2] = -0.880501;

    M_inletRadiusRef = 1.219238;
    M_outletRadiusRef1 = 0.504580;
    M_outletRadiusRef2 = 0.555306;

    resetInletOutlets();
}

void
Aorta::
resetInletOutlets()
{
    GeometricFace inlet(M_inletCenterRef, M_inletNormalRef, M_inletRadiusRef, 1, -1);
    GeometricFace outlet1(M_outletCenterRef1, M_outletNormalRef1, M_outletRadiusRef1, 2, -1);
    GeometricFace outlet2(M_outletCenterRef2, M_outletNormalRef2, M_outletRadiusRef2, 3, -1);

    M_inlet = inlet;
    M_outlets.clear();
    M_outlets.push_back(outlet1);
    M_outlets.push_back(outlet2);
}


std::string
Aorta::
getOptionalParameter(unsigned int index)
{
}

void
Aorta::
applyNonAffineTransformation(bool transformMesh)
{
}

}
