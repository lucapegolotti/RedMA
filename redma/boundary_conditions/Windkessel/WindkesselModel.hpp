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

#ifndef WINDKESSELMODEL_HPP
#define WINDKESSELMODEL_HPP

#include <redma/boundary_conditions/aBCModel.hpp>
#include <redma/boundary_conditions/Windkessel/WindkesselPressureDrop.hpp>

namespace RedMA
{

// For a reference see:
// "The nested block preconditioning technique for the incompressible Navier-Stokes
// equations with emphasis on haemodynamic simulations" - Liu, Yang, Dong, Marsden
class WindkesselModel : public aBCModel
{
public:
    WindkesselModel(const DataContainer& data, const std::string& dataEntry,
                    const unsigned int& indexOutlet);

    virtual double getNeumannCondition(const double& time, const double& rate) override;

private:
    double                                              M_C;  // compliance
    double                                              M_Rp; // proximal resistance
    double                                              M_Rd; // distal resistance
};

}

#endif // WINDKESSELMODEL_HPP
