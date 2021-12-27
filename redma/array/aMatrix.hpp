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

#ifndef aMATRIX_HPP
#define aMATRIX_HPP

#include <redma/RedMA.hpp>
#include <redma/array/aDataWrapper.hpp>
#include <redma/array/aVector.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/array/TypesUtils.hpp>

#include <memory>

namespace RedMA
{

/// Abstract class for a matrix.
class aMatrix : public aDataWrapper
{
public:

    /// Empty constructor.
    aMatrix() : M_nRows(0), M_nCols(0) {}

    /// Destructor.
    virtual ~aMatrix() {};

    /*! \brief Addition operator.
     *
     * \param other The matrix to add.
     */
    virtual void add(shp<aMatrix> other) = 0;

    /*! \brief Multiplication by scalar operator.
     *
     * \param other The coefficient to multiply.
     */
    virtual void multiplyByScalar(const double& coeff) = 0;

    /*! \brief Multiplication by vector operator.
     *
     * \param other Shared pointer to the vector to multiply.
     * \return vector Shared pointer to the output vector.
     */
    virtual shp<aVector> multiplyByVector(shp<aVector> vector) = 0;

    /*! \brief Multiplication by matrix operator.
     *
     * \param other Shared pointer to the matrix to multiply.
     * \return vector Shared pointer to the output matrix.
     */
    virtual shp<aMatrix> multiplyByMatrix(shp<aMatrix> other) = 0;

    /*! \brief Transpose operator.
     *
     * \return Transposed matrix.
     */
    virtual shp<aMatrix> transpose() const = 0;

    /*! \brief Save content of the wrapper to file.
     *
     * \param namefile Name of the output file.
     */
    virtual void dump(std::string filename) const = 0;

    /*! \brief Getter for the number of rows.
     *
     * \return Number of rows.
     */
    inline unsigned int nRows() const {return M_nRows;}

    /*! \brief Getter for the number of cols.
     *
     * \return Number of cols.
     */
    inline unsigned int nCols() const {return M_nCols;}

    /*! \brief Compute norm inf of the matrix.
     *
     * \return The infinite norm.
     */
    virtual double normInf() = 0;

    /*! \brief Compute Frobenius norm of the matrix.
     *
     * \return The infinite norm.
     */
    virtual double normFrobenius() = 0;

    /*! \brief Getter for a block.
     *
     * This function only makes sense in case of block matrices, but is included
     * here to simplify the access to blocks in certain parts of the code. If not
     * overloaded, an exception is thrown.
     *
     * \param row The desired row.
     * \param col The desired col.
     * \return Shared pointer to the block.
     */
    virtual shp<aMatrix> block(const unsigned int& row,
                               const unsigned int& col) const
    {
        throw new Exception("block(row,col) function not overloaded!");
    }

protected:
    unsigned int        M_nRows;
    unsigned int        M_nCols;
};

}

#endif // aMATRIX_HPP
