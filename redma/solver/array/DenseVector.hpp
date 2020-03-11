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

#include <redma/solver/array/aMatrix.hpp>
#include <redma/utils/Exception.hpp>

#include <Epetra_SerialDenseVector.h>

#include <fstream>

#define DENSEVECTOR         Epetra_SerialDenseVector

namespace RedMA
{

class DenseVector : public aMatrix
{
public:
    DenseVector();

    DenseVector operator+(const DenseVector& other);

    DenseVector operator-(const DenseVector& other);

    DenseVector& operator+=(const DenseVector& other);

    DenseVector& operator-=(const DenseVector& other);

    DenseVector& operator*=(const double& coeff);

    DenseVector& operator=(const std::shared_ptr<DENSEVECTOR>& other);

    void hardCopy(const DenseVector& other);

    void softCopy(const DenseVector& other);

    double norm2() const;

    std::shared_ptr<DENSEVECTOR>& data();

    std::shared_ptr<DENSEVECTOR> data() const;

    std::string getString(const char& delimiter) const;

    unsigned int getNumRows() const;

    void dump(std::string filename) const;

private:
    std::shared_ptr<DENSEVECTOR>  M_vector;
};

}

#endif // DENSEVECTOR_HPP
