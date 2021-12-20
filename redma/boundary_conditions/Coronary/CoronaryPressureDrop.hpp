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
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#ifndef CORONARYPRESSUREDROP_HPP
#define CORONARYPRESSUREDROP_HPP

#include <redma/boundary_conditions/aPressureDrop.hpp>

namespace RedMA
{

class CoronaryPressureDrop : public aPressureDrop {

public:
    CoronaryPressureDrop(const double& Ca, const double& Cim,
                         const double& Ram, const double& Rvm,
                         const double& Rv);

    virtual shp<aVector> getZeroVector() const override;

    virtual shp<aMatrix> getMass(const double& time,
                                 const shp<aVector>& sol) override;

    virtual shp<aMatrix> getMassJacobian(const double& time,
                                         const shp<aVector>& sol) override;

    virtual shp<aVector> getRightHandSide(const double& time,
                                          const shp<aVector>& sol) override;

    virtual shp<aMatrix> getJacobianRightHandSide(const double& time,
                                                  const shp<aVector>& sol) override;

    void setIntramyocardialPressure(const double& time,
                                    std::function<double(double)> Pim) {M_Pim = Pim(time);};

private:
    double                                M_Ca;
    double                                M_Cim;
    double                                M_Ram;
    double                                M_Rvm;
    double                                M_Rv;
    double                                M_Pim;

};

} // namespace RedMA


#endif //CORONARYPRESSUREDROP_HPP
