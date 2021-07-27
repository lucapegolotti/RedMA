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

#ifndef BLOCKVECTOR_HPP
#define BLOCKVECTOR_HPP

#include <redma/utils/Exception.hpp>
#include <redma/array/DenseVector.hpp>
#include <redma/array/DistributedVector.hpp>
#include <redma/array/Wrap.hpp>

#include <boost/numeric/ublas/matrix.hpp>

#include <math.h>

namespace RedMA
{

/// Class handling a block vector.
class BlockVector : public aVector
{
    typedef boost::numeric::ublas::matrix<shp<aVector>>        Grid;

public:

    /// Default constructor.
    BlockVector();

    /// Default destructor.
    virtual ~BlockVector() {};

    /*! \brief Constructor.
     *
     * \param nRows Number of rows.
     */
    BlockVector(const unsigned int& nRows);

    /*! \brief Addition operator.
     *
     * The compatibility of the input is checked internally.
     *
     * \param other The vector to add.
     */
    virtual void add(shp<aVector> other) override;

    /*! \brief Multiplication by scalar operator.
     *
     * \param coeff The coefficient to multiply.
     */
    virtual void multiplyByScalar(const double& coeff) override;

    /*! \brief Save content of the BlockVector to file.
     *
     * This method is not implemented and raises an exception.
     *
     * \param namefile Name of the output file.
     */
    virtual void dump(std::string namefile) const override;

    /*! \brief Shallow copy.
     *
     * The compatibility of the input is checked internally.
     * The shared pointers to the data of the argument are copied.
     *
     * \param other Shared pointer to another BlockVector.
     */
    virtual void shallowCopy(shp<aDataWrapper> other) override;

    /*! \brief Deep copy.
     *
     * The compatibility of the input is checked internally.
     * The data of the argument are copied.
     *
     * \param other Shared pointer to another BlockVector.
     */
    virtual void deepCopy(shp<aDataWrapper> other) override;

    /*! \brief Clones the block vector.
     *
     * \return Raw pointer to a copy of the block vector.
     */
    virtual BlockVector* clone() const override;

    /*! \brief Returns true if all the blocks are zero.
     *
     * \return True if all the blocks are zero.
     */
    virtual bool isZero() const override;

    /*! \brief Returns a string with all the components of the vector.
     *
     * This method is not implemented and raises an exception.
     *
     * \param delimiter The delimiter used to separate the components.
     * \return The desired string.
     */
    std::string getString(const char& delimiter) const override;

    /*! \brief Norm 2 of the vector.
     *
     * \return The norm.
     */
    double norm2() const override;

    /*! \brief Setter for a block.
     *
     * \param row The desired row.
     * \param vectpr Shared pointer to the vector to set (shallow).
     */
    void setBlock(const unsigned int& iblock,
                  shp<aVector> vector);

    /*! \brief Check if all the vectors are of a given type.
     *
     * \param The type.
     * \return True if all the blocks are of the given type.
     */
    bool globalTypeIs(Datatype type);

    /*! \brief Maximum depth of the block vectpr.
     *
     * If the vectpr is composed only of nonblock vectprs (e.g., Dense vectors),
     * the level is 1, if it comprises at least one block vector with level 1,
     * it is of level 2, and so on.
     *
     * \return The level.
     */
    unsigned int level();

    /*! \brief Converts internal blocks to a given type.
     *
     * Returns a block vector in which all the internal blocks are converted to
     * a given type. At the moment, we only support the conversion from
     * dense to sparse.
     *
     * \param type The Datatype.
     * \param comm The MPI Communicator.
     * \return The block vector with converted data types.
     */
    shp<BlockVector> convertInnerTo(Datatype type,
                                    shp<Epetra_Comm> comm = nullptr);

    /*! \brief Getter for a block.
     *
     * \param iblock The desired row.
     * \return Shared pointer to the block.
     */
    shp<aVector> block(const unsigned int& iblock) const override;

    /*! \brief Get a subvector.
     *
     * The blocks are not copied.
     *
     * \param ibegin Row index of the top block.
     * \param iend Row index of the bottom block.
     * \return Shared pointer to the subvector.
     */
    shp<BlockVector> getSubvector(const unsigned int& ibegin,
                                  const unsigned int& iend) const;

    /*! \brief Change dimensions of the block vector.
     *
     * \param nRows Number of rows.
     */
    void resize(const unsigned int& nRows);

    /// Set the internal MPI Communicator by looking at the internal blocks.
    void findComm();

    /*! \brief Getter for the MPI Communicator.
     *
     * \return The MPI Communicator.
     */
    shp<Epetra_Comm> commPtr() {findComm(); return M_comm;}

    /*! \brief Access operator.
     *
     * Method not implemented. It raises an exception.
     *
     * \param The index of the component to access.
     * \return The value of the desired component.
     */
    virtual double operator()(unsigned int index) override
    {
        throw new Exception("operator() undefined for BlockVector");
    }

    /*! \brief Returns a nullptr.
     *
     * This method does not makes sense for block matrices.
     *
     * \return A nullptr.
     */
    virtual shp<void> data() const override {return nullptr;};

    /*! \brief Set the internal data.
     *
     * This method throws an exception as it does not make sense for block matrices.
     *
     * \param data The data to be set.
     */
    virtual void setData(shp<void>) override
    {
        throw new Exception("setData undefined for BlockVector");
    };

    /*! \brief Getter for the type.
     *
     * \return Returns BLOCK.
     */
    virtual Datatype type() const override {return BLOCK;}

protected:
    shp<Epetra_Comm>              M_comm;
    Grid                          M_vectorGrid;
};

}

#endif // BLOCKVECTOR_HPP
