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

#ifndef BLOCKMATRIX_HPP
#define BLOCKMATRIX_HPP

#include <Exception.hpp>

#include <AbstractMatrix.hpp>

#include <Matrix.hpp>
#include <BlockVector.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapVector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace RedMA
{

// eventually this will become abstract with respect to type of inner Matrices
class BlockMatrix : public AbstractMatrix
{
    typedef Epetra_CrsMatrix                                     MatrixCrs;
    typedef std::shared_ptr<Epetra_CrsMatrix>                    MatrixCrsPtr;
    typedef boost::numeric::ublas::matrix<MatrixCrsPtr >         MatrixCrsGrid;
    typedef LifeV::MatrixEpetra<double>                          MatrixEpetra;
    typedef std::shared_ptr<MatrixEpetra>                        MatrixEpetraPtr;
    typedef boost::numeric::ublas::matrix<Matrix<MatrixEpetra> > MatrixGrid;
    typedef LifeV::VectorEpetra                                  VectorEpetra;
    typedef std::shared_ptr<VectorEpetra>                        VectorEpetraPtr;
    typedef LifeV::MapEpetra                                     Map;
    typedef std::shared_ptr<Map>                                 MapPtr;

public:
    BlockMatrix();

    BlockMatrix(unsigned int numRows, unsigned int numCols);

    // attention: I am switching this copy constructor to soft copy but this may
    // potentially break everything
    BlockMatrix(const BlockMatrix& other);

    void resize(unsigned int numRows, unsigned int numCols);

    unsigned int getNumberRows();

    unsigned int getNumberCols();

    Matrix<MatrixEpetra>& block(unsigned int row, unsigned int col);

    Matrix<MatrixEpetra> block(unsigned int row, unsigned int col) const;

    void add(const BlockMatrix& other);

    BlockVector operator*(const BlockVector& vector);

    void copyBlock(unsigned int rows, unsigned int cols,
                   MatrixEpetraPtr matrix);

    MapPtr rangeMap(unsigned int row, unsigned int col) const;

    MapPtr domainMap(unsigned int row, unsigned int col) const;

    BlockMatrix& operator*=(const double& coeff);

    MatrixCrsGrid getGrid();

    // void spy();

    // void singleNorms1();

    // void setMaps(std::vector<MapPtr> rangeMaps, std::vector<MapPtr> domainMaps);

    void printPattern();

private:
    MatrixCrsGrid        M_grid;
    MatrixGrid           M_matrixGrid;
    unsigned int         M_rows;
    unsigned int         M_cols;
    std::vector<MapPtr>  M_rangeMaps;
    std::vector<MapPtr>  M_domainMaps;
    MapPtr               M_globalDomainMap;
    MapPtr               M_globalRangeMap;
    Matrix<MatrixEpetra> M_mmm;
};

}  // namespace RedMA

#endif  // BLOCKMATRIX_HPP
