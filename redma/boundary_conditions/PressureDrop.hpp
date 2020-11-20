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

#ifndef PRESSUREDROP_HPP
#define PRESSUREDROP_HPP

#include <redma/RedMA.hpp>

#include <redma/solver/time_marching_algorithms/aFunctionProvider.hpp>
#include <redma/array/Double.hpp>
#include <redma/array/BlockVector.hpp>
#include <redma/array/BlockMatrix.hpp>

#include <fstream>

namespace RedMA
{

// we compute the pressure drop for the three element windkessel model following
// "The nested block preconditioning technique for the incompressible Navier-Stokes
// equations with emphasis on hemodynamic simulations" - Liu, Yang, Dong, Marsden
class PressureDrop : public aFunctionProvider
{
public:

    PressureDrop(const double& C, const double& Rp, const double& Rd);

    virtual shp<aVector> getZeroVector() const override;

    virtual shp<aMatrix> getMass(const double& time,
                                const shp<aVector>& sol) override;

    virtual shp<aMatrix> getMassJacobian(const double& time,
                                        const shp<aVector>& sol) override;

    virtual shp<aVector> getRightHandSide(const double& time,
                                         const shp<aVector>& sol) override;

    virtual shp<aMatrix> getJacobianRightHandSide(const double& time,
                                                 const shp<aVector>& sol) override;

    inline void setFlowRate(const double& Q) {M_Q = Q;}

    virtual void apply0DirichletBCs(shp<aVector> vector) const override {}

    virtual void applyDirichletBCs(const double& time, shp<aVector> vector) const override {}

    void setExtrapolatedSolution(const shp<aVector>& exSol) override {throw new Exception("function must still be implemented PressureDrop");}

private:
    double                      M_C;  // compliance
    double                      M_Rp; // proximal resistance
    double                      M_Rd; // distal resistance
    double                      M_Q;  // flow rate
};

}

#endif // PRESSUREDROP_HPP
