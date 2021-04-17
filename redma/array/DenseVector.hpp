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

#ifndef DENSEVECTOR_HPP
#define DENSEVECTOR_HPP

#include <redma/array/aVector.hpp>
#include <redma/utils/Exception.hpp>

#include <Epetra_SerialDenseVector.h>
#include <lifev/core/array/VectorEpetra.hpp>

#include <fstream>
#include <memory>

namespace RedMA
{

/// Wrapper of a dense vector.
class DenseVector : public aVector
{
public:

    /// Default constructor.
    DenseVector();

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

    /*! \brief Save content of the DenseVector to file.
     *
     * \param namefile Name of the output file.
     */
    virtual void dump(std::string namefile) const override;

    /*! \brief Shallow copy.
     *
     * The compatibility of the input is checked internally.
     * The shared pointer to the data of the argument is copied.
     *
     * \param other Shared pointer to another DenseVector.
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

    /*! \brief Clones the dense vector.
     *
     * \return Raw pointer to a copy of the dense vector.
     */
    virtual DenseVector* clone() const override;

    /*! \brief Returns true if the internal vector is not set.
     *
     * \return True if the internal vector is not set.
     */
    virtual bool isZero() const override;

    /*! \brief Norm 2 of the vector.
     *
     * \return The norm.
     */
    double norm2() const override;

    /*! \brief Returns a string with all the components of the vector.
     *
     * \param delimiter The delimiter used to separate the components.
     * \return The desired string.
     */
    std::string getString(const char& delimiter) const override;

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

    /*! \brief Setter for the internal vector.
     *
     * \param matrix Shared pointer to the vector.
     */
    void setVector(shp<DENSEVECTOR> vector);

    /*! \brief Getter for the internal vector.
     *
     * \return Shared pointer to the internal vector.
     */
    shp<DENSEVECTOR> getVector();

    /*! \brief Convert the dense vector to epetra.
     *
     * \param comm The MPI Communicator
     * \return Shared pointer to the epetra vector
     */
    shp<VECTOREPETRA> toVectorEpetraPtr(shp<Epetra_Comm> comm) const;

    /*! \brief Access operator.
     *
     * \param The index of the component to access.
     * \return The value of the desired component.
     */
    virtual double operator()(unsigned int index) override {return M_vector->operator()(index);}

    /*! \brief Getter for the type.
     *
     * \return Returns DENSE.
     */
    virtual Datatype type() const override {return DENSE;}

private:
    shp<DENSEVECTOR>        M_vector;
};

}

#endif // DENSEVECTOR_HPP
