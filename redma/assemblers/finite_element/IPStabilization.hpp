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

#ifndef IPSTABILIZATION_HPP
#define IPSTABILIZATION_HPP

#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/VectorElemental.hpp>
#include <lifev/core/fem/AssemblyElemental.hpp>
#include <lifev/core/fem/Assembly.hpp>

#include <redma/assemblers/finite_element/NavierStokesStabilization.hpp>

namespace RedMA {

/*! \brief Class for the stabilization of the Navier-Stokes equations with IP method.
 */
class IPStabilization : public NavierStokesStabilization
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
    IPStabilization(const DataContainer& data,
                    shp<FESPACE> fespaceVelocity,
                    shp<FESPACE> fespacePressure,
                    shp<ETFESPACE3> etfespaceVelocity,
                    shp<ETFESPACE1> etfespacePressure,
                    EPETRACOMM comm);

    /*! \brief Setup method.
     *
     */
    virtual void setup() override;

    /*! \brief Assemble and get the mass.
     *
     * \param sol The solution.
     * \param rhs The right hand side.
     * \return The desired matrix.
     */
    virtual shp<BlockMatrix> getMass(shp<BlockVector> sol = nullptr,
                                     shp<BlockVector> rhs = nullptr) override;

    /*! \brief Assemble and get the mass Jacobian.
     *
     * \param sol The solution.
     * \param rhs The right hand side.
     * \return The desired matrix.
     */
    virtual shp<BlockMatrix> getMassJacobian(shp<BlockVector> sol = nullptr,
                                             shp<BlockVector> rhs = nullptr) override;

    /*! \brief Assemble and get the Jacobian.
    *
    * \param sol The solution.
    * \param rhs The right hand side.
    * \return The desired matrix.
    */
    virtual shp<BlockMatrix> getJacobian(shp<BlockVector> sol = nullptr,
                                         shp<BlockVector> rhs = nullptr) override;

    /*! \brief Assemble and get the residual.
     *
     * \param sol The solution.
     * \param rhs The right hand side.
     * \return The desired vector.
     */
    virtual shp<BlockVector> getResidual(shp<BlockVector> sol,
                                         shp<BlockVector> rhs = nullptr) override;

private:

    shp<MATRIXEPETRA> stabilizePressure();

    double                               M_gamma_pressure;

    shp<LifeV::CurrentFE>                M_pressureFESpaceSide1;
    shp<LifeV::CurrentFE>                M_pressureFESpaceSide2;

};

} // namespace RedMA

#endif //IPSTABILIZATION_HPP
