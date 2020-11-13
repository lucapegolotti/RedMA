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

#define DENSEVECTOR         Epetra_SerialDenseVector

namespace RedMA
{

class DenseVector : public aVector
{
public:
    DenseVector();

    DenseVector(const DenseVector& vector);

    virtual void add(std::shared_ptr<aVector> other) override;

    virtual void multiplyByScalar(const double& coeff) override;

    virtual void dump(std::string namefile) const override;

    virtual void softCopy(std::shared_ptr<aVector> other) override;

    virtual void hardCopy(std::shared_ptr<aVector> other) override;

    virtual aVector* cloneVector() const override;

    virtual bool isZero() override;

    double norm2() const override;

    std::string getString(const char& delimiter) const override;

    virtual std::shared_ptr<void> data() const override;

    virtual void setData(std::shared_ptr<void> data) override;

    void setVector(std::shared_ptr<DENSEVECTOR> vector);

    std::shared_ptr<LifeV::VectorEpetra> toVectorEpetraPtr(std::shared_ptr<Epetra_Comm> comm) const;

    virtual double operator()(unsigned int index) override {return M_vector->operator()(index);}

private:
    std::shared_ptr<DENSEVECTOR>  M_vector;
};

}

#endif // DENSEVECTOR_HPP
