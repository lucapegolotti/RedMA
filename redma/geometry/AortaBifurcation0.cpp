#include "AortaBifurcation0.hpp"

namespace RedMA
{

AortaBifurcation0::
AortaBifurcation0(commPtr_Type comm, std::string name, bool verbose) :
  BuildingBlock(comm, "normal", verbose)
{
    M_name = name;
    M_datafileName = "data_mesh";
    M_meshName = "others/aortabif0_fine.mesh";

    // center of inlet (reference configuration)
    M_inletCenterRef[0] = -0.0340;
    M_inletCenterRef[1] = 3.0347;
    M_inletCenterRef[2] = 11.3465;

    // center of outlet (reference configuration)
    M_outletCenterRef1[0] = -5.716039;
    M_outletCenterRef1[1] = 5.621785;
    M_outletCenterRef1[2] = -20.529217;

    // center of outlet (reference configuration)
    M_outletCenterRef2[0] = 8.266322;
    M_outletCenterRef2[1] = 4.993420;
    M_outletCenterRef2[2] = -19.989317;

    // center of outlet (reference configuration)
    M_outletCenterRef3[0] = -5.716039;
    M_outletCenterRef3[1] = 5.621785;
    M_outletCenterRef3[2] = -20.529217;

    // center of outlet (reference configuration)
    M_outletCenterRef4[0] = 8.266322;
    M_outletCenterRef4[1] = 4.993420;
    M_outletCenterRef4[2] = -19.989317;

    // normal of inlet (reference configuration)
    M_inletNormalRef[0] = -0.128;
    M_inletNormalRef[1] = -0.9803;
    M_inletNormalRef[2] = 0.1504257624;

    // outlet of outlet (reference configuration)
    M_outletNormalRef1[0] = -0.318503;
    M_outletNormalRef1[1] = 0.325307;
    M_outletNormalRef1[2] = -0.890355;

    // outlet of outlet (reference configuration)
    M_outletNormalRef2[0] = 0.272467;
    M_outletNormalRef2[1] = 0.387918;
    M_outletNormalRef2[2] = -0.880501;

    // outlet of outlet (reference configuration)
    M_outletNormalRef3[0] = -0.318503;
    M_outletNormalRef3[1] = 0.325307;
    M_outletNormalRef3[2] = -0.890355;

    // outlet of outlet (reference configuration)
    M_outletNormalRef4[0] = 0.272467;
    M_outletNormalRef4[1] = 0.387918;
    M_outletNormalRef4[2] = -0.880501;

    M_inletRadiusRef = 1.570332;
    M_outletRadiusRef1 = 0.504580;
    M_outletRadiusRef2 = 0.555306;
    M_outletRadiusRef3 = 0.504580;
    M_outletRadiusRef4 = 0.555306;

    M_wallFlag = 1;
    resetInletOutlets();
}

void
AortaBifurcation0::
resetInletOutlets()
{
    GeometricFace inlet(M_inletCenterRef, M_inletNormalRef, M_inletRadiusRef, 3, -1);
    GeometricFace outlet1(M_outletCenterRef1, M_outletNormalRef1, M_outletRadiusRef1, 2, -1);
    GeometricFace outlet2(M_outletCenterRef2, M_outletNormalRef2, M_outletRadiusRef2, 4, -1);
    GeometricFace outlet3(M_outletCenterRef3, M_outletNormalRef3, M_outletRadiusRef3, 5, -1);
    GeometricFace outlet4(M_outletCenterRef4, M_outletNormalRef4, M_outletRadiusRef4, 6, -1);

    M_inlet = inlet;
    M_outlets.clear();
    M_outlets.push_back(outlet1);
    M_outlets.push_back(outlet2);
    M_outlets.push_back(outlet3);
    M_outlets.push_back(outlet4);
}


std::string
AortaBifurcation0::
getOptionalParameter(unsigned int index)
{
}

void
AortaBifurcation0::
applyNonAffineTransformation(bool transformMesh)
{
}

}
