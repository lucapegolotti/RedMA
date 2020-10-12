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

#include <lifev/core/array/VectorEpetra.hpp>

#define VECTOREPETRA        LifeV::VectorEpetra

namespace RedMA
{

class DistributedVector : public aVector
{
public:
    DistributedVector();

    DistributedVector(const DistributedVector& vector);

    virtual void add(std::shared_ptr<aVector> other) override;

    virtual void multiplyByScalar(const double& coeff) override;

    virtual void dump(std::string namefile) const override;

    virtual void softCopy(std::shared_ptr<aVector> other) override;

    virtual void hardCopy(std::shared_ptr<aVector> other) override;

    virtual aVector* clone() const override;

    virtual bool isZero() override;

    std::string getString(const char& delimiter) const override;

    double norm2() const override;

    virtual std::shared_ptr<void> data() const override;

    virtual void setData(std::shared_ptr<void> data) override;

    double maxMagnitude3D() const;

    DenseVector toDenseVector() const;

    std::shared_ptr<Epetra_Comm> commPtr() {return M_vector->mapPtr()->commPtr();}

    std::shared_ptr<DenseVector> toDenseVectorPtr() const;

    static std::shared_ptr<DistributedVector> convertDenseVector(
        std::shared_ptr<DenseVector> denseVector,
        std::shared_ptr<Epetra_Comm> comm);

    void setVector(std::shared_ptr<VECTOREPETRA> vector);

private:
    std::shared_ptr<VECTOREPETRA>  M_vector;
};

std::shared_ptr<DistributedVector> epetraToDistributed(std::shared_ptr<VECTOREPETRA> vector);

}

#endif // DISTRIBUTEDVECTOR_HPP
