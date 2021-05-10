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

#ifndef BASISFUNCTIONFUNCTOR_HPP
#define BASISFUNCTIONFUNCTOR_HPP

#include <functional>

#include <lifev/core/array/VectorSmall.hpp>

#include <redma/geometry/building_blocks/BuildingBlock.hpp>

namespace RedMA
{

/*! \brief A functor implementing a Lagrange multiplier basis function.
 *
 * These basis functions are used to enforce the coupling across interfaces.
 * The basis functions are indexed, such that when the evaluation operators are
 * called the returned function values corresponds to the current function
 * index.
 */
class BasisFunctionFunctor
{
public:
    typedef double                                            return_Type;
    typedef LifeV::VectorSmall<3>                             Vector3D;
    typedef std::function<double(double const&,
                                 double const&,
                                 double const&,
                                 double const&,
                                 unsigned int const& )>       Function;

    /*! \brief Constructor taking a GeometricFace.
     *
     * \param face The GeometricFace where the basis functions are defined.
     */
    BasisFunctionFunctor(const GeometricFace& face);

    /*! \brief Evaluation operator.
     *
     * \param pos Position where the function has to be evaluated.
     * \return Value of the basis function.
     */
    virtual return_Type operator()(const Vector3D& pos) = 0;

    /*! \brief Create a std::function out of evaluateOperator.
     *
     * The function is of type std::function<double(double const&,
     * double const&, double const&, double const&, unsigned int const& )>.
     * \return The function.
     */
    Function function();

    /*! \brief Evaluate the operator at given coordinates.
     *
     * The t and index variables are actually unused and only included for
     * compatibility with LifeV.
     * \param t Unused variable.
     * \param x x coordinate.
     * \param y y coordinate.
     * \param z z coordinate.
     * \param index Unused variable.
     */
    return_Type evaluateOperator(const double& t,
                                 const double& x,
                                 const double& y,
                                 const double& z,
                                 unsigned int const& index);

    /* \brief Set the index of the current basis function.
     *
     * \param index The index.
     */
    virtual void setIndex(const unsigned int& index);

    /* \brief Get the number of basis functions.
     *
     * \return The number of basis functions
     */
    unsigned int getNumBasisFunctions() const;

    /* \brief Get the type of basis function as a string.
     *
     * \return The type of the basis function.
     */
    std::string getType() const {return M_type;};

protected:
    void getLocalXAndY(const Vector3D& pos,
                       double& x,
                       double& y);

    void getThetaAndRadius(const Vector3D& pos,
                           double& theta,
                           double& radius);

    GeometricFace           M_face;
    unsigned int            M_index;
    unsigned int            M_nBasisFunctions;
    // versor to compute the normal
    Vector3D                M_e;
    Vector3D                M_eOrth;
    std::string             M_type;
};

}  // namespace RedMA

#endif  // BASISFUNCTIONFUNCTOR_HPP
