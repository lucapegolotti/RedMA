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

#ifndef CLOTHFUNCTIONFUNCTOR_HPP
#define CLOTHFUNCTIONFUNCTOR_HPP


#include <functional>

#include <lifev/core/array/VectorSmall.hpp>

#include <redma/geometry/building_blocks/BuildingBlock.hpp>

namespace RedMA
{

/*! \brief A functor implementing the indicator function for the cloth modeling.
*
*/
class ClothFunctionFunctor
{
public:
    typedef double                                            return_Type;
    typedef LifeV::VectorSmall<3>                             Vector3D;
    typedef LifeV::MatrixSmall<3,3>                           Matrix3D;
    typedef std::function<double(double const&,
            double const&,
            double const&,
            double const&,
            unsigned int const& )>       Function;

    /*! \brief Constructor taking the cloth center and radius
     *
     * \param center Center of the cloth
     * \param radius Radius of the cloth
     */
    ClothFunctionFunctor(const LifeV::Vector3D& center, const double& radius,
                         const LifeV::Vector3D& normal, const LifeV::Vector3D& tangent,
                         const LifeV::Vector3D& shape_coefficients);

    /*! \brief Evaluation operator.
     *
     * \param pos Position where the function has to be evaluated.
     * \return Value of the cloth indicator function.
     */
    return_Type operator()(const Vector3D& pos);

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

protected:

    /*! \brief Evaluation operator.
     *
     * \param pos Position where the function has to be evaluated.
     * \return Value of the cloth indicator function.
     */
    return_Type norm(const Vector3D& pos) const;

    Vector3D                  M_center;
    double                    M_radius;
    Vector3D                  M_normal;
    Vector3D                  M_tangent1;
    Vector3D                  M_tangent2;
    Vector3D                  M_shapeCoeffs;
    Matrix3D                  M_normMatrix;
};

}  // namespace RedMA


#endif //CLOTHFUNCTIONFUNCTOR_HPP
