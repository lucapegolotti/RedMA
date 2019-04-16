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
    M_parametersMap["r0"] = M_inletRadiusRef;
    M_parametersMap["r1"] = M_outletRadiusRef;
}

}
