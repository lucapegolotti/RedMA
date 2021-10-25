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

    double tol = 1e-10;
    M_isResistanceBC = (std::abs(M_C) < tol) && (std::abs(M_Rd) < tol);
    if (!M_isResistanceBC)
    {
        M_pressureDrop.reset(new WindkesselPressureDrop(M_C, M_Rd));

        M_bdf.reset(new BDF(data, M_pressureDrop));
        M_bdf->setLinearSolver();
    }
    else
        printlog(GREEN, "[WindkesselModel] The Windkessel model reduces to a Resistance  "
                        "one as the capacitance and the distal resistance have been set to 0.");
}

double
WindkesselModel::
getNeumannCondition(const double& time, const double& rate)
{
    double retVal;
    if (!M_isResistanceBC)
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

        retVal = ramp * M_Rp * rate +
                 spcast<DoubleVector>(M_pressureDropSolution->block(0))->getValue() +
                 M_Pd(time);
    }
    else
    {
        retVal = M_Rp * rate +
                 M_Pd(time);
    }

    return retVal;
}

double
WindkesselModel::
getNeumannJacobian(const double &time, const double &rate)
{
    if (M_isResistanceBC)
        return M_Rp;
    else
        return aBCModel::getNeumannJacobian(time, rate);

}

void
WindkesselModel::
shiftSolutions()
{
    if (!M_isResistanceBC)
        aBCModel::shiftSolutions();
}

}
