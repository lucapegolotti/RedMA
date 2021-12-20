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

#ifndef DOUBLE_HPP
#define DOUBLE_HPP

#include <iostream>
#include <redma/array/aVector.hpp>

namespace RedMA
{

/// Wrapper of a double.
class DoubleVector : public aVector
{
public:

    /// Default constructor.
    DoubleVector();

    /*! \brief Addition operator.
     *
     * The compatibility of the input is checked internally.
     *
     * \param other The double to add.
     *//*
    virtual void add(shp<aMatrix> other) override;*/

    /*! \brief Multiplication by scalar operator.
     *
     * \param coeff The coefficient to multiply.
     */
    virtual void multiplyByScalar(const double& coeff) override;

    /*! \brief Multiplication by another DoubleVector.
     *
     * The compatibility of the input is checked internally.
     *
     * \param other Shared pointer to the DoubleVector to multiply.
     * \return vector Shared pointer to the output.
     *//*
    virtual shp<aVector> multiplyByVector(shp<aVector> vector) override;*/

    /*! \brief Multiplication by another DoubleVector.
     *
     * The compatibility of the input is checked internally.
     *
     * \param other Shared pointer to the DoubleVector to multiply.
     * \return vector Shared pointer to the output.
     *//*
    virtual shp<aMatrix> multiplyByMatrix(shp<aMatrix> other) override;*/

    /// This method is not implemented.
    virtual void dump(std::string namefile) const override;

    /*! \brief Return true if the internal double is zero.
     *
     * \return True if the internal double is zero.
     */
    virtual bool isZero() const override;

    /// This method is not implemented.
    virtual DoubleVector* clone() const override;

    /*! \brief Returns the value of the double regardless of the input.
     *
     * \param index This parameter is not used.
     */
    virtual double operator()(unsigned int index) override;

    /*! \brief Addition operator.
     *
     * The compatibility of the input is checked internally.
     *
     * \param other The double to add.
     */
    virtual void add(shp<aVector> other) override;

    /*! \brief Shallow copy.
     *
     * The compatibility of the input is checked internally.
     * The shared pointer to the data of the argument is copied.
     *
     * \param other Shared pointer to another DoubleVector.
     */
    virtual void shallowCopy(shp<aDataWrapper> other) override;

    /*! \brief Deep copy.
     *
     * The compatibility of the input is checked internally.
     * The data of the argument is copied.
     *
     * \param other Shared pointer to another DoubleVector.
     */
    virtual void deepCopy(shp<aDataWrapper> other) override;

    /// This method is not implemented.
    virtual std::string getString(const char& delimiter) const override;

    /*! \brief Returns the absolute value of the double.
     *
     * \return The absolute value of the double.
     */
    virtual double norm2() const override {return std::abs(M_double);}

    /*! \brief Setter for the internal double.
     *
     * \param data The double.
     */
    void setValue(double data) {M_double = data;}

    /*! \brief Getter for the internal double.
     *
     * \return The internal double.
     */
    double getValue() {return M_double;}

    /*! \brief Getter for the type.
     *
     * \return Returns DOUBLE.
     */
    virtual Datatype type() const override {return DOUBLE;}

    /*! \brief Returns a shared pointer containing the double.
     *
     * \return The shared pointer to the double.
     */
    virtual shp<void> data() const override;

    /// Method not implemented.
    virtual void setData(shp<void>) override {};

    /// Method not implemented.
    // virtual shp<aMatrix> transpose() const override {return nullptr;};

    /*! \brief Compute norm inf of the double.
     *
     * \return The infinite norm.
     *//*
    virtual double normInf() override {return std::abs(M_double);}*/

private:
    double  M_double;
};

// specification of template function to avoid ambiguous cast
template<>
shp<DoubleVector> convert(shp<aDataWrapper> container);

}

#endif // DOUBLE_HPP
