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
#include <lifev/eta/expression/InteriorIntegrate.hpp>

#include <functional>

namespace RedMA
{

/*! \brief Abstract class for the stabilization of the Navier-Stokes equations.
 */
class NavierStokesStabilization
{
public:
    /*! \brief Constructor.
     *
     * \param data A DataContainer object.
     * \param fespaceVelocity Finite element space for the velocity.
     * \param fespacePressure Finite element space for the pressure.
     * \param etfespaceVelocity ETA finite element space for the velocity.
     * \param etfespacePressure ETA finite element space for the pressure.
     */
    NavierStokesStabilization(const DataContainer& data,
                              shp<FESPACE> fespaceVelocity,
                              shp<FESPACE> fespacePressure,
                              shp<ETFESPACE3> etfespaceVelocity,
                              shp<ETFESPACE1> etfespacePressure,
                              EPETRACOMM comm);

    /*! \brief Setup method.
     *
     */
    virtual void setup() = 0;

    /*! \brief Setter for the density and the velocity.
     *
     * \param density The value of the density.
     * \param viscosity The value of the viscosity.
     */
    void setDensityAndViscosity(const double& density,
                                const double& viscosity);

    /*! \brief Assemble and get the mass.
     *
     * \param sol The solution.
     * \param rhs The right hand side.
     * \return The desired matrix.
     */
    virtual shp<BlockMatrix> getMass(shp<BlockVector> sol,
                                     shp<BlockVector> rhs) = 0;

    /*! \brief Assemble and get the mass Jacobian.
     *
     * \param sol The solution.
     * \param rhs The right hand side.
     * \return The desired matrix.
     */
    virtual shp<BlockMatrix> getMassJacobian(shp<BlockVector> sol,
                                             shp<BlockVector> rhs) = 0;

    /*! \brief Assemble and get the Jacobian.
     *
     * \param sol The solution.
     * \param rhs The right hand side.
     * \return The desired matrix.
     */
    virtual shp<BlockMatrix> getJacobian(shp<BlockVector> sol,
                                         shp<BlockVector> rhs) = 0;

    /*! \brief Assemble and get the residual.
     *
     * \param sol The solution.
     * \param rhs The right hand side.
     * \return The desired vector.
     */
    virtual shp<BlockVector> getResidual(shp<BlockVector> sol,
                                         shp<BlockVector> rhs) = 0;

protected:
    EPETRACOMM                          M_comm;
    bool                                M_verbose;

    double                              M_density;
    double                              M_viscosity;
    shp<FESPACE>                        M_velocityFESpace;
    shp<FESPACE>                        M_pressureFESpace;
    shp<ETFESPACE3>                     M_velocityFESpaceETA;
    shp<ETFESPACE1>                     M_pressureFESpaceETA;
    shp<BlockMatrix>                    M_jac;
    shp<BlockMatrix>                    M_massJac;
    shp<BlockMatrix>                    M_mass;

    std::string                         M_velocityOrder;
};

}  // namespace RedMA

#endif  // NAVIERSTOKESSTABILIZATION_HPP
