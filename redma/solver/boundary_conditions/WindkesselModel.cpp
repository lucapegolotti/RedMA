#include "WindkesselModel.hpp"

namespace RedMA
{

WindkesselModel::
WindkesselModel(const DataContainer& data, const std::string& dataEntry,
                const unsigned int& indexOutlet)
{
    M_C = data(dataEntry + "/C", 0);
    M_Rp = data(dataEntry + "/Rp", 0);
    M_Rd = data(dataEntry + "/Rd", 0);
    // attention: here we assume that the timestep is constant!
    M_dt = data("time_discretization/dt", 0.01);
    M_pressureDrop.reset(new PressureDrop(M_C, M_Rp, M_Rd));
    M_Pd = data.getDistalPressure(indexOutlet);
    CoutRedirecter ct;
    ct.redirect();
    M_bdf.reset(new BDF<Double, Double>(data, M_pressureDrop));
    ct.restore();
}

double
WindkesselModel::
getNeumannCondition(const double& time, const double& rate)
{
    M_pressureDrop->setFlowRate(rate);
    M_pressureDropSolution.resize(1);

    int status = -1;

    CoutRedirecter ct;
    ct.redirect();
    M_pressureDropSolution = M_bdf->advance(time, M_dt, status);
    ct.restore();

    if (status)
        throw new Exception("Error in WindkesselModel: status != 0");

    BlockVector<Double> sol(1);
    sol.block(0).data() = M_Rp * rate + M_pressureDropSolution.block(0).data() + M_Pd(time);

    return sol.block(0).data();
}

double
WindkesselModel::
getNeumannJacobian(const double& rate)
{
    return M_Rp + 1.0/M_C * M_dt * M_bdf->getRhsCoeff();
}

void
WindkesselModel::
shiftSolutions()
{
    M_bdf->shiftSolutions(M_pressureDropSolution);
}

}
