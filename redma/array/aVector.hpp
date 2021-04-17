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

#ifndef aVECTOR_HPP
#define aVECTOR_HPP

#include <redma/RedMA.hpp>
#include <redma/array/aDataWrapper.hpp>
#include <redma/array/Datatypes.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/array/TypesUtils.hpp>

#include <memory>

namespace RedMA
{

/// Abstract class for a vector.
class aVector : public aDataWrapper
{
public:
    /// Empty constructor.
    aVector() : M_nRows(0) {}

    /// Destructor.
    virtual ~aVector() {};

    /*! \brief Addition operator.
     *
     * \param other The vector to add.
     */
    virtual void add(shp<aVector> other) = 0;

    /*! \brief Multiplication by scalar operator.
     *
     * \param other The coefficient to multiply.
     */
    virtual void multiplyByScalar(const double& coeff) = 0;

    /*! \brief Returns a string with all the components of the vector.
     *
     * \param delimiter The delimiter used to separate the components.
     * \return The desired string.
     */
    virtual std::string getString(const char& delimiter) const = 0;

    /*! \brief Norm 2 of the vector.
     *
     * \return The norm.
     */
    virtual double norm2() const = 0;

    /*! \brief Access operator.
     *
     * \param The index of the component to access.
     * \return The value of the desired component.
     */
    virtual double operator()(unsigned int index) = 0;

    /*! \brief Getter for the number of rows.
     *
     * \return Number of rows.
     */
    inline unsigned int nRows() const {return M_nRows;}

    /*! \brief Getter for a block.
     *
     * This function only makes sense in case of block vectors, but is included
     * here to simplify the access to blocks in certain parts of the code. If not
     * overloaded, an exception is thrown.
     *
     * \param row The desired row.
     * \return Shared pointer to the block.
     */
    virtual shp<aVector> block(const unsigned int& row) const
    {
        throw new Exception("block(row) function not overloaded!");
    }

protected:
    unsigned int        M_nRows;
};

}

#endif // aVECTOR_HPP
