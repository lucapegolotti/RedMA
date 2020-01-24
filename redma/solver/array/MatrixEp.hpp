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

#ifndef MATRIXEP_HPP
#define MATRIXEP_HPP

#include <redma/solver/array/aMatrix.hpp>
#include <redma/solver/array/VectorEp.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>

#define MATRIXEPETRA        LifeV::MatrixEpetra<double>

namespace RedMA
{

class MatrixEp : public aMatrix
{
public:
    MatrixEp();

    MatrixEp operator+(const MatrixEp& other);

    MatrixEp& operator+=(const MatrixEp& other);

    MatrixEp& operator*=(const double& coeff);

    VectorEp operator*(const VectorEp& vector);

    void hardCopy(const MatrixEp& other);

    void softCopy(const MatrixEp& other);

    void getRowProperty(std::shared_ptr<LifeV::MapEpetra>& outMap);

    void getColProperty(std::shared_ptr<LifeV::MapEpetra>& outMap);

    std::shared_ptr<MATRIXEPETRA>& data();

private:
    std::shared_ptr<MATRIXEPETRA>  M_matrix;
};

}

#endif // MATRIXEP_HPP
