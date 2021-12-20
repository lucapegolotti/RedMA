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

#ifndef DENSEMATRIX_HPP
#define DENSEMATRIX_HPP

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <redma/utils/PrintLog.hpp>
#include <redma/array/aMatrix.hpp>
#include <redma/array/DenseVector.hpp>

#include <Epetra_SerialDenseMatrix.h>

#include <fstream>

namespace RedMA
{

/// Wrapper of a dense matrix.
class DenseMatrix : public aMatrix
{
public:

    /// Default constructor.
    DenseMatrix();

    /*! \brief Addition operator.
     *
     * The compatibility of the input is checked internally.
     *
     * \param other The matrix to add.
     */
    virtual void add(shp<aMatrix> other) override;

    /*! \brief Multiplication by scalar operator.
     *
     * \param coeff The coefficient to multiply.
     */
    virtual void multiplyByScalar(const double& coeff) override;

    /*! \brief Multiplication by matrix operator.
     *
     * The compatibility of the input is checked internally.
     *
     * \param other Shared pointer to the matrix to multiply.
     * \return vector Shared pointer to the output matrix.
     */
    virtual shp<aMatrix> multiplyByMatrix(shp<aMatrix> other) override;

    /*! \brief Transpose operator.
     *
     * \return Transposed matrix.
     */
    virtual shp<aMatrix> transpose() const override;

    /*! \brief Multiplication by vector operator.
     *
     * The compatibility of the input is checked internally.
     *
     * \param other Shared pointer to the vector to multiply.
     * \return vector Shared pointer to the output vector.
     */
    virtual shp<aVector> multiplyByVector(shp<aVector> vector) override;

    /*! \brief Shallow copy.
     *
     * The compatibility of the input is checked internally.
     * The shared pointer to the data of the argument is copied.
     *
     * \param other Shared pointer to another DenseMatrix.
     */
    virtual void shallowCopy(shp<aDataWrapper> other) override;

    /*! \brief Deep copy.
     *
     * The compatibility of the input is checked internally.
     * The data of the argument is copied.
     *
     * \param other Shared pointer to another DenseMatrix.
     */
    virtual void deepCopy(shp<aDataWrapper> other) override;

    /*! \brief Returns true if the internal matrix is not set.
     *
     * \return True the internal matrix is not set.
     */
    virtual bool isZero() const override;

    /*! \brief Getter for the data.
     *
     * \return Shared pointer to the data.
     */
    virtual shp<void> data() const override;

    /*! \brief Setter for the data.
     *
     * \param data Shared pointer to the data.
     */
    virtual void setData(shp<void> data) override;

    /*! \brief Save content of the DenseMatrix to file.
     *
     * \param namefile Name of the output file.
     */
    virtual void dump(std::string filename) const override;

    /*! \brief Clones the dense matrix.
     *
     * \return Raw pointer to a copy of the dense matrix.
     */
    virtual DenseMatrix* clone() const override;

    /*! \brief Getter for the type.
     *
     * \return Returns DENSE.
     */
    virtual Datatype type() const override {return DENSE;}

    /*! \brief Setter for the internal matrix.
     *
     * \param matrix Shared pointer to the matrix.
     */
    void setMatrix(shp<DENSEMATRIX> matrix);

    /*! \brief Getter for the internal matrix.
     *
     * \return Shared pointer to the internal matrix.
     */
    shp<DENSEMATRIX> getMatrix();

    /*! \brief Compute norm inf of the matrix.
     *
     * \return The infinite norm.
     */
    virtual double normInf() override {return M_matrix->NormInf();};

    /*! \brief Compute Frobenius of the matrix.
     *
     * \return The Frobenius norm.
     */
    virtual double normFrobenius() override {return 0.0;};  // not implemented for dense matrices!

private:

    shp<DENSEMATRIX>        M_matrix;
};

}

#endif // DENSEMATRIX_HPP
