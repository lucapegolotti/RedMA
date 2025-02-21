#include "Aorta.hpp"

namespace RedMA
{

Aorta::
Aorta(EPETRACOMM comm, std::string refinement, bool verbose, bool add_rings) :
  BuildingBlock(comm, refinement, verbose),
  M_add_rings(add_rings)
{
    M_name = "aorta";
    M_datafileName = "data_mesh";

    if (!std::strcmp(refinement.c_str(), "coarse"))
        M_meshName = add_rings ? "others/aorta_coarse_rings.mesh" : "others/aorta_coarse.mesh";
    else if (!std::strcmp(refinement.c_str(), "normal"))
        M_meshName = add_rings ? "others/aorta_normal_rings.mesh" : "others/aorta_normal.mesh";
    else if (!std::strcmp(refinement.c_str(), "fine"))
        M_meshName = add_rings ? "others/aorta_rings.mesh" :  "others/aorta.mesh";
    else
        throw new Exception("Undefined refinement: " + refinement);

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

    //M_inletRadiusRef = 1.219238;
    M_outletRadiusRef1 = 0.504580;
    M_outletRadiusRef2 = 0.555306;

    M_inletRadius1Ref = 1.354801;
    M_inletRadius2Ref = 1.072587;
    M_inletTangent1Ref[0] = -0.05504388;
    M_inletTangent1Ref[1] = 0.91404855;
    M_inletTangent1Ref[2] = 0.40185248;
    M_inletTangent2Ref[0] = -0.89172613;
    M_inletTangent2Ref[1] = 0.13606545;
    M_inletTangent2Ref[2] = -0.43163724;

    M_wallFlag = 10;

    resetInletOutlets();
}

void
Aorta::
resetInletOutlets()
{
    // elliptic inlet face
    GeometricFace inlet(M_inletCenterRef, M_inletNormalRef,
                        M_inletRadius1Ref, M_inletRadius2Ref, M_inletTangent1Ref, M_inletTangent2Ref,
                        1, M_add_rings ? 1000 : -1);
    GeometricFace outlet1(M_outletCenterRef1, M_outletNormalRef1, M_outletRadiusRef1, 2,
                          M_add_rings ? 2000 : -1);
    GeometricFace outlet2(M_outletCenterRef2, M_outletNormalRef2, M_outletRadiusRef2, 3,
                          M_add_rings ? 3000 : -1);

    M_inlets.clear();
    M_inlets.push_back(inlet);
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
