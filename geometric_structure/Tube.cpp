#include "Tube.hpp"

namespace ReMA
{

Tube::Tube(commPtr_Type comm, bool verbose) :
  BuildingBlock(comm, verbose)
{
    M_name = "Tube";

    M_datafileName = "tube_coarse_data";

    // it is important to fill parametersMap right at this level because then
    // the keys will be used in the parser to check the values in the XML file
    // center of inlet
    M_parametersMap["x0"] = 0.0;
    M_parametersMap["y0"] = 0.0;
    M_parametersMap["z0"] = 0.0;

    // center of outlet
    M_parametersMap["x1"] = 0.0;
    M_parametersMap["y1"] = 0.0;
    M_parametersMap["z1"] = 15.0;

    // normal of inlet
    M_parametersMap["i0"] = 0.0;
    M_parametersMap["j0"] = 0.0;
    M_parametersMap["k0"] = -1.0;

    // normal of outlet
    M_parametersMap["i1"] = 0.0;
    M_parametersMap["j1"] = 0.0;
    M_parametersMap["k1"] = 1.0;

    // radius of inlet
    M_parametersMap["r0"] = 3.0;

    // radius of outlet
    M_parametersMap["r1"] = 3.0;
}

}
