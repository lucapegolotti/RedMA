#include "AortaBifurcation1.hpp"

namespace RedMA
{

AortaBifurcation1::
AortaBifurcation1(commPtr_Type comm, std::string name, bool verbose) :
  BuildingBlock(comm, "normal", verbose)
{
    M_name = name;
    M_datafileName = "data_mesh";
    M_meshName = "others/aortabif1_fine.mesh";

    // center of inlet (reference configuration)
    M_inletCenterRef[0] = 1.389511;
    M_inletCenterRef[1] = 1.291508;
    M_inletCenterRef[2] = 13.856705;

    M_inletNormalRef[0] = -0.509002;
    M_inletNormalRef[1] = -0.010880;
    M_inletNormalRef[2] = -0.860697;

    // center of outlet (reference configuration)
    M_outletCenterRef1[0] = 1.568659;
    M_outletCenterRef1[1] = 1.587641;
    M_outletCenterRef1[2] = 15.215374;

    M_outletNormalRef1[0] = 0.049329;
    M_outletNormalRef1[1] = 0.054728;
    M_outletNormalRef1[2] = 0.997282;

    // center of outlet (reference configuration)
    M_outletCenterRef2[0] = 2.636050;
    M_outletCenterRef2[1] = 1.221660;
    M_outletCenterRef2[2] = 15.60280;

    // center of outlet (reference configuration)
    M_outletNormalRef2[0] = 0.715379;
    M_outletNormalRef2[1] = -0.190554;
    M_outletNormalRef2[2] = 0.672252;

    M_inletRadiusRef = 0.441261;
    M_outletRadiusRef1 = 0.347327;
    M_outletCenterRef2 = 0.344803392;

    M_wallFlag = 1;
    resetInletOutlets();
}

void
AortaBifurcation1::
resetInletOutlets()
{
    GeometricFace inlet(M_inletCenterRef, M_inletNormalRef, M_inletRadiusRef, 2, -1);
    GeometricFace outlet1(M_outletCenterRef1, M_outletNormalRef1, M_outletRadiusRef1, 3, -1);
    GeometricFace outlet2(M_outletCenterRef2, M_outletNormalRef2, M_outletRadiusRef2, 4, -1);

    M_inlet = inlet;
    M_outlets.clear();
    M_outlets.push_back(outlet1);
    M_outlets.push_back(outlet2);
}


std::string
AortaBifurcation1::
getOptionalParameter(unsigned int index)
{
}

void
AortaBifurcation1::
applyNonAffineTransformation(bool transformMesh)
{
}

}
