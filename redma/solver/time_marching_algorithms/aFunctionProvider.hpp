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

#ifndef aFUNCTIONPROVIDER_HPP
#define aFUNCTIONPROVIDER_HPP

#include <redma/RedMA.hpp>

#include <redma/array/aVector.hpp>
#include <redma/array/aMatrix.hpp>

#include <fstream>

namespace RedMA
{

/*! \brief Generic function provider.
 *
 * We are interested in solving equations in the form
 *  \f[
 *    M(u) \dot{u} = F(t,u).
 *  \f]
 * The classes inheriting from this abstract class need to implement all the
 * elements necessary to the solution of the equation with the Newton method.
 */
class aFunctionProvider
{
public:
    /// Empty constructor.
    aFunctionProvider() {};

    /*! \brief Getter for the zero vector.
     *
     * \return Shared pointer to the zero vector.
     */
    virtual shp<aVector> getZeroVector() const = 0;

    /*! \brief Getter for the mass matrix.
     *
     * \param time The current time.
     * \param Shared pointer to the solution.
     * \return Shared pointer to the mass matrix.
     */
    virtual shp<aMatrix> getMass(const double& time,
                                 const shp<aVector>& sol) = 0;

    /*! \brief Getter for the pressure mass matrix.
     *
     * \param time The current time.
     * \param Shared pointer to the solution.
     * \return Shared pointer to the pressure mass matrix.
     */
    virtual shp<aMatrix> getPressureMass(const double& time, 
                                         const shp<aVector>& sol) = 0;

    /*! \brief Getter for the mass matrix Jacobian.
     *
     * \param time The current time.
     * \param Shared pointer to the solution.
     * \return Shared pointer to the mass matrix Jacobian.
     */
    virtual shp<aMatrix> getMassJacobian(const double& time,
                                         const shp<aVector>& sol) = 0;

    /*! \brief Getter for the right-hand side.
     *
     * \param time The current time.
     * \param Shared pointer to the solution.
     * \return Shared pointer to the right-hand side.
     */
    virtual shp<aVector> getRightHandSide(const double& time,
                                          const shp<aVector>& sol) = 0;

    /*! \brief Getter for the right-hand side Jacobian.
     *
     * \param time The current time.
     * \param Shared pointer to the solution.
     * \return Shared pointer to the right-hand side Jacobian.
     */
    virtual shp<aMatrix> getJacobianRightHandSide(const double& time,
                                                  const shp<aVector>& sol) = 0;

    /*! \brief Apply homogeneous Dirichlet boundary conditions to a vector.
     *
     * \param vector Shared pointer to the vector.
     */
    virtual void apply0DirichletBCs(shp<aVector> vector) const = 0;

    /*! \brief Apply Dirichlet boundary conditions to a vector.
     *
     * \param time The current time.
     * \param vector Shared pointer to the vector.
     */
    virtual void applyDirichletBCs(const double& time,
                                   shp<aVector> vector) const = 0;

    /*! \brief Set the extrapolated solution.
     *
     * \param exSol Shared pointer to the extrapolated solution.
     */
    virtual void setExtrapolatedSolution(const shp<aVector>& exSol) = 0;
};

}

#endif // aFUNCTIONPROVIDER_HPP
