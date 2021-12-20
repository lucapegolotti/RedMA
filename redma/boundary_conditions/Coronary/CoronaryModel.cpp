#include "CoronaryModel.hpp"

namespace RedMA {


CoronaryModel::
CoronaryModel(const DataContainer& data, const std::string& dataEntry,
              const unsigned int& indexOutlet) :
    aBCModel(data, dataEntry, indexOutlet)
{
    M_Ca = data(dataEntry + "/Ca", 0.0);
    M_Cim = data(dataEntry + "/Cim", 0.0);
    M_Ra = data(dataEntry + "/Ra", 0.0);
    M_Ram = data(dataEntry + "/Ram", 0.0);
    M_Rvm = data(dataEntry + "/Rvm", 0.0);
    M_Rv = data(dataEntry + "/Rv", 0.0);

    M_pressureDrop.reset(new CoronaryPressureDrop(M_Ca, M_Cim, M_Ram, M_Rvm, M_Rv));

    M_bdf.reset(new BDF(data, M_pressureDrop));
    M_bdf->setLinearSolver();

    M_Pim = data.getIntramyocardialPressure();
}

double
CoronaryModel::
getNeumannCondition(const double& time, const double& rate)
{
    M_pressureDrop->setFlowRate(rate);
    spcast<CoronaryPressureDrop>(M_pressureDrop)->setIntramyocardialPressure(time - M_dt, M_Pim);

    int status = -1;

    // time - M_dt because we are assuming that we are already in a Newton iteration
    M_pressureDropSolution = spcast<BlockVector>(M_bdf->advance(time - M_dt, M_dt, status));

    if (status)
        throw new Exception("Error in Coronary Model: status != 0");

    // we do this to smooth out the value of Ra, which otherwise at the beginning
    // of the simulation is too high
    double ramp = M_data.evaluateRamp(time);

    double retVal = ramp * M_Ra * rate +
                    spcast<DoubleVector>(M_pressureDropSolution->block(0))->getValue() +
                    M_Pd(time);

    return retVal;
}

}  // namespace RedMA