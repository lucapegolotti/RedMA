#include "WindkesselModel.hpp"

namespace RedMA
{

WindkesselModel::
WindkesselModel(const DataContainer& data, const std::string& dataEntry,
                const unsigned int& indexOutlet)
{
    M_C = data(dataEntry + "/C", 0.0);
    M_Rp = data(dataEntry + "/Rp", 0.0);
    M_Rd = data(dataEntry + "/Rd", 0.0);
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
    // time - M_dt because we are assuming that we are already in a newton iteration
    M_pressureDropSolution = M_bdf->advance(time - M_dt, M_dt, status);
    ct.restore();

    if (status)
        throw new Exception("Error in WindkesselModel: status != 0");

    return M_Rp * rate + M_pressureDropSolution.block(0).data() * 0 + M_Pd(time);
}

double
WindkesselModel::
getNeumannJacobian(const double& time, const double& rate)
{
    const double epsabs = 1e-8;
    const double epsrel = 1e-5;

    double eps = epsabs > epsrel * std::abs(rate) ? epsabs : epsrel * std::abs(rate);

    std::cout << "eps = " << eps << std::endl << std::flush;

    // we approximate the jacobian via finite differences
    double jac1 = getNeumannCondition(time, rate + eps / 2.);
    double jac2 = getNeumannCondition(time, rate - eps / 2.);
    // this is to restore M_pressureDropSolution (otherwise shiftSolutions will
    // utilize a wrong one)
    getNeumannCondition(time, rate);
    std::cout << (jac1 - jac2) / eps << std::endl << std::flush;
    return (jac1 - jac2) / eps;
}

void
WindkesselModel::
shiftSolutions()
{
    M_bdf->shiftSolutions(M_pressureDropSolution);
}

}
