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

#ifndef BLOCKMATRIX_HPP
#define BLOCKMATRIX_HPP

#include <redma/utils/Exception.hpp>
#include <redma/array/aMatrix.hpp>
#include <redma/array/BlockVector.hpp>
#include <redma/array/DistributedVector.hpp>
#include <redma/array/DenseVector.hpp>
#include <redma/array/SparseMatrix.hpp>
#include <redma/array/DenseMatrix.hpp>
#include <redma/array/DoubleMatrix.hpp>
#include <redma/array/Wrap.hpp>

#include <boost/numeric/ublas/matrix.hpp>

#include <lifev/core/array/MapEpetra.hpp>


namespace RedMA
{

/// Class handling a block matrix.
class BlockMatrix : public aMatrix
{

    typedef boost::numeric::ublas::matrix<shp<aMatrix>>       Grid;

public:

    /// Default constructor.
    BlockMatrix();

    /// Default destructor.
    virtual ~BlockMatrix() {};

    /*! \brief Constructor.
     *
     * \param nRows Number of rows.
     * \param nCols Number of cols.
     */
    BlockMatrix(const unsigned int& nRows,
                const unsigned int& nCols);

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
     * The shared pointers to the data of the argument are copied.
     *
     * \param other Shared pointer to another BlockMatrix.
     */
    virtual void shallowCopy(shp<aDataWrapper> other) override;

    /*! \brief Deep copy.
     *
     * The compatibility of the input is checked internally.
     * The data of the argument are copied.
     *
     * \param other Shared pointer to another BlockMatrix.
     */
    virtual void deepCopy(shp<aDataWrapper> other) override;

    /*! \brief Returns true if all the blocks are zero.
     *
     * \return True if all the blocks are zero.
     */
    virtual bool isZero() const override;

    /*! \brief Save content of the BlockMatrix to file.
     *
     * At the end of the output filename, we append _block_ and the indices of
     * the block.
     *
     * \param namefile Name of the output file.
     */
    virtual void dump(std::string filename) const override;

    /*! \brief Clones the block matrix.
     *
     * \return Raw pointer to a copy of the block matrix.
     */
    virtual BlockMatrix* clone() const override;

    /*! \brief Change dimensions of the block matrix.
     *
     * \param nRows Number of rows.
     * \param nCols Number of cols.
     */
    void resize(const unsigned int& nRows,
                const unsigned int& nCols);

    /*! \brief Getter for a block.
     *
     * \param iblock The desired row.
     * \param jblock The desired col.
     * \return Shared pointer to the block.
     */
    shp<aMatrix> block(const unsigned int& iblock,
                       const unsigned int& jblock) const override;

    /*! \brief Setter for a block.
     *
     * \param iblock The desired row.
     * \param jblock The desired col.
     * \param matrix Shared pointer to the matrix to set (shallow).
     */
    void setBlock(const unsigned int& iblock,
                  const unsigned int& jblock,
                  shp<aMatrix> matrix);

    /*! \brief Get a submatrix.
     *
     * The blocks are not copied.
     *
     * \param ibegin Row index of the top left block.
     * \param iend Row index of the bottom right block.
     * \param jbegin Column index of the top left block.
     * \param jend Column index of the bottom right block.
     * \return Shared pointer to the submatrix.
     */
    shp<BlockMatrix> getSubmatrix(const unsigned int& ibegin,
                                  const unsigned int& iend,
                                  const unsigned int& jbegin,
                                  const unsigned int& jend) const;

    /*! \brief Maximum depth of the block matrix.
     *
     * If the matrix is composed only of nonblock matrices (e.g., Dense matrices),
     * the level is 1, if it comprises at least one block matrix with level 1,
     * it is of level 2, and so on.
     *
     * \return The level.
     */
    unsigned int level();

    /*! \brief Check if all the matrices are of a given type.
     *
     * \param The type.
     * \return True if all the blocks are of the given type.
     */
    bool globalTypeIs(Datatype type);

    /*! \brief Converts internal blocks to a given type.
     *
     * Returns a block matrix in which all the internal blocks are converted to
     * a given type. At the moment, we only support the conversion from
     * dense to sparse.
     *
     * \param type The Datatype.
     * \param comm The MPI Communicator.
     * \return The block matrix with converted data types.
     */
    shp<BlockMatrix> convertInnerTo(Datatype type,
                                    shp<Epetra_Comm> comm = nullptr);

    /// Print pattern of the matrx.
    void printPattern() const;

    /// Set the internal MPI Communicator by looking at the internal blocks.
    void findComm();

    /*! \brief Getter for the MPI Communicator.
     *
     * \return The MPI Communicator.
     */
    shp<Epetra_Comm> commPtr() {findComm(); return M_comm;}

    /*! \brief Getter for the type.
     *
     * \return Returns BLOCK.
     */
    virtual Datatype type() const override {return BLOCK;}

    virtual shp<void> data() const override {return nullptr;};

    /*! \brief Set the internal data.
     *
     * This method throws an exception as it does not make sense for block matrices.
     *
     * \param data The data to be set.
     */
    virtual void setData(shp<void> data) override
    {
        throw new Exception("setData undefined for BlockMatrix");
    };

    /*! \brief Compute norm inf of the matrix.
     *
     * \return The infinite norm.
     */
    virtual double normInf() override;

    /*! \brief Compute teh Frobenius norm of the matrix.
     *
     * \return The Frobenius norm.
     */
    virtual double normFrobenius() override;

protected:

    shp<Epetra_Comm>              M_comm;
    Grid                          M_matrixGrid;
};

shp<SparseMatrix> blockMatrixToSparseMatrix(shp<BlockMatrix> matrix);

}

#endif // BLOCKMATRIX_HPP
