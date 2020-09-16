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

#ifndef NAVIERSTOKESSTABILIZATION_HPP
#define NAVIERSTOKESSTABILIZATION_HPP

// this is taken from StabilizationSUPG_semi_implicit
// MACRO TO DEFINE TAU_M
#define TAU_M 	       value(1.0)/( eval(squareroot,TAU_M_DEN) )
#define TAU_M_DEN      TAU_M_DEN_DT + TAU_M_DEN_VEL + TAU_M_DEN_VISC
#define TAU_M_DEN_DT   value(M_density*M_density)*value(M_timeOrder*M_timeOrder)/value(M_dt * M_dt)
#define TAU_M_DEN_VEL  value(M_density*M_density)*dot(value(M_velocityFESpaceETA, *velocityRep), G*value(M_velocityFESpaceETA, *velocityRep))
#define TAU_M_DEN_VISC value(M_C_I)*value(M_viscosity*M_viscosity)*dot(G,G)

#define TAU_C ( value(1.0)/( dot(g, TAU_M * g ) ) )

#define VH M_velocityFESpaceETA,*velocityRep
#define PH M_pressureFESpaceETA,*pressureRep
#define MOMENTUM_R value(M_density) * value(VH) * grad(VH) - value(M_density) * value(M_velocityFESpaceETA,*velocityRhsRep) + grad(PH) - value(M_viscosity) * laplacian(VH)
#define MOMENTUM_R_DER value(M_density) * phi_j * grad(VH) +  value(M_density) * value(VH) * grad(phi_j) - value(M_viscosity) * laplacian(phi_j)

#include <redma/RedMA.hpp>

#include <redma/array/BlockVector.hpp>
#include <redma/array/BlockMatrix.hpp>
#include <redma/array/DistributedVector.hpp>
#include <redma/array/SparseMatrix.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/problem/DataContainer.hpp>

#include <lifev/eta/utils/Functions.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <functional>

namespace RedMA
{

class NavierStokesStabilization
{
public:
    NavierStokesStabilization(const DataContainer& data,
                              SHP(FESPACE) fespaceVelocity,
                              SHP(FESPACE) fespacePressure,
                              SHP(ETFESPACE3) etfespaceVelocity,
                              SHP(ETFESPACE1) etfespacePressure);

    void setDensityAndViscosity(const double& density, const double& viscosity);

    virtual BlockMatrix getMass(const BlockVector& sol,
                                const BlockVector& rhs) = 0;

    virtual BlockMatrix getMassJac(const BlockVector& sol,
                                   const BlockVector& rhs) = 0;

    virtual BlockMatrix getJac(const BlockVector& sol,
                                         const BlockVector& rhs) = 0;

    virtual BlockVector getResidual(const BlockVector& sol,
                                              const BlockVector& rhs) = 0;

protected:
    double                          M_density;
    double                          M_viscosity;
    unsigned int                    M_timeOrder;
    unsigned int                    M_velocityOrder;
    SHP(FESPACE)                    M_velocityFESpace;
    SHP(FESPACE)                    M_pressureFESpace;
    SHP(ETFESPACE3)                 M_velocityFESpaceETA;
    SHP(ETFESPACE1)                 M_pressureFESpaceETA;
    SHP(BlockMatrix)                M_jac;
    SHP(BlockMatrix)                M_massJac;
    SHP(BlockMatrix)                M_mass;
    double                          M_dt;
    double                          M_C_I;
};

}  // namespace RedMA

#endif  // NAVIERSTOKESSTABILIZATION_HPP
