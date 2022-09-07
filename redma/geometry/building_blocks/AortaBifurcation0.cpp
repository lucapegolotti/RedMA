#include "AortaBifurcation0.hpp"

namespace RedMA
{

AortaBifurcation0::
AortaBifurcation0(EPETRACOMM comm, std::string refinement,
                  std::string name, bool verbose) :
  BuildingBlock(comm, "normal", verbose)
{
    M_name = name;
    M_datafileName = "data_mesh";
    M_meshName = "others/aortabif0_" + refinement + ".mesh";

    // center of inlet (reference configuration)
    M_inletCenterRef[0] = -0.074727;
    M_inletCenterRef[1] = 2.963938;
    M_inletCenterRef[2] = 11.383856;

    // center of outlet (reference configuration)
    M_outletCenterRef1[0] = -0.653389;
    M_outletCenterRef1[1] = -0.420303;
    M_outletCenterRef1[2] = 10.520454;

    // center of outlet (reference configuration)
    M_outletCenterRef2[0] = -1.001712;
    M_outletCenterRef2[1] = 0.132182;
    M_outletCenterRef2[2] = 11.836387;

    // center of outlet (reference configuration)
    M_outletCenterRef3[0] = -0.497380;
    M_outletCenterRef3[1] = 0.955569;
    M_outletCenterRef3[2] = 12.644601;

    // center of outlet (reference configuration)
    M_outletCenterRef4[0] = 0.587853;
    M_outletCenterRef4[1] = 1.496366;
    M_outletCenterRef4[2] = 12.628265;

    // normal of inlet (reference configuration)
    M_inletNormalRef[0] = -0.088701;
    M_inletNormalRef[1] = 0.989459;
    M_inletNormalRef[2] = -0.114466;

    // normal of outlet (reference configuration)
    M_outletNormalRef1[0] = -0.158101;
    M_outletNormalRef1[1] = -0.817230;
    M_outletNormalRef1[2] = -0.554201;

    // normal of outlet (reference configuration)
    M_outletNormalRef2[0] = -0.418575;
    M_outletNormalRef2[1] = -0.286510;
    M_outletNormalRef2[2] = 0.861804;

    // normal of outlet (reference configuration)
    M_outletNormalRef3[0] = -0.358876;
    M_outletNormalRef3[1] = -0.285907;
    M_outletNormalRef3[2] = 0.888519;

    // normal of outlet (reference configuration)
    M_outletNormalRef4[0] = 0.560367;
    M_outletNormalRef4[1] = -0.284151;
    M_outletNormalRef4[2] = 0.777976;

    M_inletRadiusRef = 1.090612;
    M_outletRadiusRef1 = 0.760686;
    M_outletRadiusRef2 = 0.371728;
    M_outletRadiusRef3 = 0.386731;
    M_outletRadiusRef4 = 0.617990;

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
    GeometricFace outlet3(M_outletCenterRef3, M_outletNormalRef3, M_outletRadiusRef3, 6, -1);
    GeometricFace outlet4(M_outletCenterRef4, M_outletNormalRef4, M_outletRadiusRef4, 5, -1);

    M_inlets.clear();
    M_inlets.push_back(inlet);
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
