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

#ifndef SPARSEMATRIX_HPP
#define SPARSEMATRIX_HPP

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <redma/utils/PrintLog.hpp>
#include <redma/array/aMatrix.hpp>
#include <redma/array/DistributedVector.hpp>
#include <redma/array/DenseMatrix.hpp>

#include <lifev/core/array/MapEpetra.hpp>

namespace RedMA
{

/// Wrapper of a sparse matrix.
class SparseMatrix : public aMatrix
{
public:

    /// Default constructor.
    SparseMatrix();

    /*! \brief Constructor.
     *
     * Build a sparse matrix out of a vector of distributed vectors.
     *
     * \param Vector of distributed vectors.
     */
    SparseMatrix(std::vector<shp<DistributedVector>> columnVectors);

    /*! \brief Constructor.
     *
     * Build a sparse matrix out of a vector of epetra vectors.
     *
     * \param Vector of vector epetra.
     */
    SparseMatrix(const std::vector<shp<VECTOREPETRA>>& columnVectors);

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
     * \param other Shared pointer to another SparseMatrix.
     */
    virtual void shallowCopy(shp<aDataWrapper> other) override;

    /*! \brief Deep copy.
     *
     * The compatibility of the input is checked internally.
     * The data of the argument is copied.
     *
     * \param other Shared pointer to another SparseMatrix.
     */
    virtual void deepCopy(shp<aDataWrapper> other) override;

    /*! \brief Returns true if the internal matrix is not set.
     *
     * \return True the internal matrix is not set.
     */
    virtual bool isZero() const override;

    /*! \brief Clones the dense matrix.
     *
     * \return Raw pointer to a copy of the sparse matrix.
     */
    virtual SparseMatrix* clone() const override;

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

    /*! \brief Save content of the SparseMatrix to file.
     *
     * \param namefile Name of the output file.
     */
    virtual void dump(std::string namefile) const override;

    /*! \brief Getter for the type.
     *
     * \return Returns SPARSE.
     */
    virtual Datatype type() const override {return SPARSE;}

    /*! \brief Convert a DenseMatrix to a SparseMatrix.
     *
     *
     * \param denseMatrix Shared pointer to the DenseMatrix.
     * \param comm The MPI Communicator.
     * \return Shared pointer to the SparseMatrix.
     */
    static shp<SparseMatrix> convertDenseMatrix(shp<DenseMatrix> denseMatrix,
                                                shp<Epetra_Comm> comm);

    /*! \brief Convert to dense matrix.
     *
     * \return The DenseMatrix.
     */
    DenseMatrix toDenseMatrix() const;

    /*! \brief Convert to a shared pointer of DenseMatrix type.
     *
     * \return The shared pointer to DenseMatrix.
     */
    shp<DenseMatrix> toDenseMatrixPtr() const;

    /*! \brief Set internal matrix.
     *
     * \param matrix shared pointer to an Epetra Matrix.
     */
    void setMatrix(shp<MATRIXEPETRA> matrix);

    /*! \brief Get the internal matrix.
     *
     * \return Shared pointer to the internal Epetra Matrix.
     */
    shp<MATRIXEPETRA> getMatrix();

    /*! \brief Getter for the MPI Communicator.
     *
     * \return The MPI Communicator.
     */
    shp<Epetra_Comm> commPtr() {return M_matrix->rangeMapPtr()->commPtr();}

    /*! \brief Compute norm inf of the matrix.
     *
     * \return The infinite norm.
     */
    virtual double normInf() override {return M_matrix->normInf();};

    /*! \brief Compute Frobenius norm of the matrix.
     *
     * \return The infinite norm.
     */
    virtual double normFrobenius() override {return M_matrix->normFrobenius();};

private:
    shp<MATRIXEPETRA>           M_matrix;
};

}

#endif // SPARSEMATRIX_HPP
