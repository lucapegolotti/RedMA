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
