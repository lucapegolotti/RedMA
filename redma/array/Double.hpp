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
#include <redma/array/aMatrix.hpp>


namespace RedMA
{

class Double : public aVector, public aMatrix
{
public:
    Double();

    virtual void add(shp<aMatrix> other) override;

    virtual void multiplyByScalar(const double& coeff) override;

    virtual shp<aVector> multiplyByVector(shp<aVector> vector) override;

    virtual shp<aMatrix> multiplyByMatrix(shp<aMatrix> other) override;

    virtual void dump(std::string namefile) const override;

    virtual bool isZero() const override;

    virtual aDataWrapper* clone() const override;

    virtual double operator()(unsigned int index) override;

    virtual void add(shp<aVector> other) override;

    virtual void shallowCopy(shp<aDataWrapper> other) override;

    virtual void deepCopy(shp<aDataWrapper> other) override;

    virtual std::string getString(const char& delimiter) const override;

    virtual double norm2() const override;

    void setValue(double data) {M_double = data;}

    double getValue() {return M_double;}

    virtual Datatype type() const override {return DOUBLE;}

    virtual shp<void> data() const override;

    virtual void setData(shp<void>) override {};

    virtual shp<aMatrix> transpose() const override {return nullptr;};

private:
    double  M_double;
};

// specification of template function to avoid ambiguous cast
template<>
shp<Double> convert(shp<aDataWrapper> container);

}

#endif // DOUBLE_HPP
