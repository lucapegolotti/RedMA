#include "aBCModel.hpp"

namespace RedMA
{

aBCModel::
aBCModel(const DataContainer& data, const std::string& dataEntry, const unsigned int& indexOutlet): M_data(data)
{
    M_Pd = data.getDistalPressure(indexOutlet);
    // attention: here we assume that the timestep is constant!
    M_dt = data("time_discretization/dt", 0.01);
}

double
aBCModel::
getNeumannJacobian(const double &time, const double &rate)
{
    const double epsabs = 1e-8;
    const double epsrel = 1e-5;

    double eps = epsabs > epsrel * std::abs(rate) ? epsabs : epsrel * std::abs(rate);

    // we approximate the jacobian via finite differences
    double jac1 = this->getNeumannCondition(time, rate + eps / 2.0);
    double jac2 = this->getNeumannCondition(time, rate - eps / 2.0);

    // this is to restore M_pressureDropSolution (otherwise shiftSolutions will utilize a wrong one)
    this->getNeumannCondition(time, rate);

    return (jac1 - jac2) / eps;
}

void
aBCModel::
shiftSolutions()
{
    M_bdf->shiftSolutions(M_pressureDropSolution);
}


} // namespace RedMA