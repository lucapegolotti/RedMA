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

#ifndef ABCMODEL_HPP
#define ABCMODEL_HPP

#include <redma/RedMA.hpp>

#include <redma/boundary_conditions/aPressureDrop.hpp>
#include <redma/solver/time_marching_algorithms/TimeMarchingAlgorithmFactory.hpp>


namespace RedMA {

class aBCModel {

public:
    aBCModel(const DataContainer& data, const std::string& dataEntry,
             const unsigned int& indexOutlet);

    virtual double getNeumannCondition(const double& time, const double& rate) = 0;

    virtual double getNeumannJacobian(const double& time, const double& rate);

    virtual void shiftSolutions();

protected:

    DataContainer                                       M_data;
    shp<BDF>                                            M_bdf;

    shp<aPressureDrop>                                  M_pressureDrop;
    shp<BlockVector>                                    M_pressureDropSolution;

    double                                              M_dt;
    std::function<double(double)>                       M_Pd; // distal reference pressure
};

} // namespace RedMA

#endif //ABCMODEL_HPP
