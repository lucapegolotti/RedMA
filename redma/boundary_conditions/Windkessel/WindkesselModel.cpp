#include "WindkesselModel.hpp"

namespace RedMA
{

WindkesselModel::
WindkesselModel(const DataContainer& data, const std::string& dataEntry,
                const unsigned int& indexOutlet) :
    aBCModel(data, dataEntry, indexOutlet)
{
    M_C = data(dataEntry + "/C", 0.0);
    M_Rp = data(dataEntry + "/Rp", 0.0);
    M_Rd = data(dataEntry + "/Rd", 0.0);

    M_pressureDrop.reset(new WindkesselPressureDrop(M_C, M_Rd));

    M_bdf.reset(new BDF(data, M_pressureDrop));
    M_bdf->setLinearSolver();
}

double
WindkesselModel::
getNeumannCondition(const double& time, const double& rate)
{
    M_pressureDrop->setFlowRate(rate);

    int status = -1;

    // time - M_dt because we are assuming that we are already in a Newton iteration
    M_pressureDropSolution = spcast<BlockVector>(M_bdf->advance(time - M_dt, M_dt, status));

    if (status)
        throw new Exception("Error in Windkessel Model: status != 0");

    // we do this to smooth out the value of Rp, which otherwise at the beginning
    // of the simulation is too high
    double ramp = M_data.evaluateRamp(time);

    double retVal = ramp * M_Rp * rate +
                    spcast<DoubleVector>(M_pressureDropSolution->block(0))->getValue() +
                    M_Pd(time);

    return retVal;
}

}
