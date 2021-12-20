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

#ifndef CORONARYMODEL_HPP
#define CORONARYMODEL_HPP

#include <redma/boundary_conditions/aBCModel.hpp>
#include <redma/boundary_conditions/Coronary/CoronaryPressureDrop.hpp>

namespace RedMA
{

class CoronaryModel : public aBCModel {

public:
    CoronaryModel(const DataContainer& data, const std::string& dataEntry,
                  const unsigned int& indexOutlet);

    virtual double getNeumannCondition(const double& time, const double& rate) override;

private:
    double                                            M_Ca;
    double                                            M_Cim;
    double                                            M_Ra;
    double                                            M_Ram;
    double                                            M_Rvm;
    double                                            M_Rv;
    std::function<double(double)>                     M_Pim;


};

} // namespace RedMA


#endif //CORONARYMODEL_HPP
