// Reduced Modeling of Arteries (RedMA)
// Copyright (C) 2019  Luca Pegolotti
//
// RedMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RedMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef INFLOWS_HPP
#define INFLOWS_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/system_solver/FunctionFunctor.hpp>

#include <cmath>
#include<Eigen/Dense>

double inflow(const double t, const std::vector<double> params, const double T);

double inflow_periodic(const double t, const std::vector<double> params, const double T, const double Tramp);

double inflow_systolic(const double t, const std::vector<double> params, const double Tramp);

double inflow_heartbeat(const double t, const std::vector<double> params, const double Tramp);

double inflow_bypass(const double t, const std::vector<double> params, const double Tramp);

struct BSpline {
    std::vector<double> knots;
    std::vector<double> controlPoints;
    unsigned int degree;
};

double splineBasisFunction(unsigned int i, unsigned int p, const std::vector<double>& knots, double t);

double evaluateBSpline(const BSpline& spline, double t);

#endif //INFLOWS_HPP
