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

#ifndef BDF_HPP
#define BDF_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/time_marching_algorithms/aTimeMarchingAlgorithm.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/array/BlockVector.hpp>
#include <redma/solver/system_solver/FunctionFunctor.hpp>
#include <redma/solver/system_solver/SystemSolver.hpp>
// #include <redma/array/DoubleVector.hpp>

#include <memory>

namespace RedMA
{

/*! \brief Backward Differentiation Formula.
 *
 * See https://en.wikipedia.org/wiki/Backward_differentiation_formula.
 */
class BDF : public aTimeMarchingAlgorithm
{
    typedef aFunctionProvider      FunProvider;

public:
    /*! \brief Constructor.
     *
     * \param datafile The DataContainer of the problem.
     */
    BDF(const DataContainer& data);

    /*! \brief Constructor.
     *
     * \param datafile The DataContainer of the problem.
     * \param funProvider The function provider.
     */
    BDF(const DataContainer& data,
        shp<FunProvider> funProvider);

    /*! \brief Constructor.
     *
     * \param datafile The DataContainer of the problem.
     * \param zeroVector Shared pointer to a zero vector.
     */
    BDF(const DataContainer& data, 
        const shp<aVector>& zeroVector);

    /*! \brief Setup function.
     *
     * \param zeroVector Shared pointer to a zero vector.
     */
    virtual void setup(const shp<aVector>& zeroVector) override;

    /*! \brief Advance function.
     *
     * \param time The time.
     * \param dt The timestep size.
     * \param status Return code; 0 if successful.
     */
    virtual shp<aVector> advance(const double& time,
                                 double& dt,
                                 int& status) override;

    /*! \brief Simple advance function used in the Membrane Model.
     *
     * \param time The time.
     * \param sol Shared pointer to the solution in order to update the displacement.
     */
    virtual shp<aVector> simpleAdvance(const double &dt, 
                                       const shp<BlockVector> &sol) override;

    /*! \brief Shift previous solutions given the new one.
     *
     * \param sol Shared pointer to the new solution.
     */
    virtual void shiftSolutions(const shp<aVector>& sol) override;

    /*! \brief Compute derivative of a function.
     *
     * \param solnp1 Solution at time n+1.
     * \param dt Timestep size.
     * \return Shared pointer to the derivative.
     */
    virtual shp<aVector> computeDerivative(const shp<aVector>& solnp1,
                                           double& dt) override;

    /*! \brief Function to compute the extrapolated solution.
     *
     * \return Shared pointer to the extrapolated solution.
     */
    virtual shp<aVector> computeExtrapolatedSolution() override;

    /*! \brief Function to combine solutions at previous time instants.
     *
     * \return Shared pointer to the combination between the previous solutions.
     */
    shp<aVector> combineOldSolutions() override;

    /*! \brief Getter of the time marching coefficients.
     *
     * \return Vector of the BDF coefficients.
     */
    std::vector<double> getCoefficients() const override;

protected:

    void setBDFCoefficients();

    void setExtrapolationCoefficients();

    std::vector<shp<BlockVector>>            M_prevSolutions;
    std::vector<double>                      M_coefficients;
    double                                   M_rhsCoeff;
    unsigned int                             M_order;
    unsigned int                             M_extrapolationOrder;
    std::vector<double>                      M_extrapolationCoefficients;
    bool                                     M_useExtrapolation;
};

}

#endif // BDF_HPP
