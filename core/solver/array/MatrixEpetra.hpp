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

#ifndef MATRIXEPETRA_HPP
#define MATRIXEPETRA_HPP

#include <Exception.hpp>

#include <Matrix.hpp>
#include <Vector.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

namespace RedMA
{

template<>
class MatrixEpetra : public Matrix<LifeV::MatrixEpetra<double> >
{
    typedef std::shared_ptr<InMatrixType>                       InMatrixPtr;
    typedef LifeV::VectorEpetra                                 VectorEpetra;
    typedef std::shared_ptr<VectorEpetra>                       VectorEpetraPtr;
    typedef LifeV::MapEpetra                                    Map;
    typedef std::shared_ptr<Map>                                MapPtr;
    typedef std::shared_ptr<LifeV::MatrixEpetra>                InMatrixEpetraPtr;

public:
    MatrixEpetra();

    MatrixEpetra(InMatrixEpetraPtr matrix);

    virtual MatrixEpetra(const MatrixEpetra& other);

    virtual void add(const MatrixEpetra& other);

    template<>
    virtual void VectorEpetra operator*(const VectorEpetra& vector);

    MatrixEpetra& operator*=(const double& coeff);

};

template<>
MatrixEpetra::
MatrixEpetra()
{

}

template<>
MatrixEpetra::
MatrixEpetra(MatrixEpetraPtr matrix)
{
    M_matrix = matrix;
}

template<>
MatrixEpetra::
MatrixEpetra(const MatrixEpetraf)

}  // namespace RedMA

#endif  // MATRIXEPETRA_HPP
