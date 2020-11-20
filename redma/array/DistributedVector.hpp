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

#ifndef DISTRIBUTEDVECTOR_HPP
#define DISTRIBUTEDVECTOR_HPP

#include <redma/array/aVector.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/array/DenseVector.hpp>

namespace RedMA
{

class DistributedVector : public aVector
{
public:
    DistributedVector();

    DistributedVector(const DistributedVector& vector);

    virtual void add(shp<aVector> other) override;

    virtual void multiplyByScalar(const double& coeff) override;

    virtual void dump(std::string namefile) const override;

    virtual void shallowCopy(shp<aDataWrapper> other) override;

    virtual void deepCopy(shp<aDataWrapper> other) override;

    virtual DistributedVector* clone() const override;

    virtual bool isZero() const override;

    std::string getString(const char& delimiter) const override;

    double norm2() const override;

    virtual shp<void> data() const override;

    virtual void setData(shp<void> data) override;

    double maxMagnitude3D() const;

    DenseVector toDenseVector() const;

    shp<Epetra_Comm> commPtr() {return M_vector->mapPtr()->commPtr();}

    shp<DenseVector> toDenseVectorPtr() const;

    static shp<DistributedVector> convertDenseVector(
        shp<DenseVector> denseVector,
        shp<Epetra_Comm> comm);

    void setVector(shp<VECTOREPETRA> vector);

    shp<VECTOREPETRA> getVector() {return M_vector;}

    virtual double operator()(unsigned int index) override {return M_vector->operator[](index);}

    virtual Datatype type() const override {return DISTRIBUTED;}

private:
    shp<VECTOREPETRA>  M_vector;
};

shp<DistributedVector> epetraToDistributed(shp<VECTOREPETRA> vector);

}

#endif // DISTRIBUTEDVECTOR_HPP
