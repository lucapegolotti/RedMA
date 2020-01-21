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

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapVector.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>
#include <lifev/core/array/MatrixEpetraStructuredUtility.hpp>

#include <boost/numeric/ublas/matrix.hpp>

namespace RedMA
{

class BlockMatrix
{
    typedef Epetra_CrsMatrix                                    Matrix;
    typedef std::shared_ptr<Matrix>                             MatrixPtr;
    typedef LifeV::MatrixEpetra<double>                         MatrixEpetra;
    typedef std::shared_ptr<MatrixEpetra>                       MatrixEpetraPtr;
    typedef boost::numeric::ublas::matrix<MatrixPtr>            Grid;
    typedef boost::numeric::ublas::matrix<MatrixEpetraPtr>      GridEpetra;
    typedef LifeV::VectorEpetra                                 VectorEpetra;
    typedef std::shared_ptr<VectorEpetra>                       VectorEpetraPtr;
    typedef LifeV::MapEpetra                                    Map;
    typedef std::shared_ptr<Map>                                MapPtr;

public:
    BlockMatrix();

    BlockMatrix(unsigned int numRows, unsigned int numCols);

    BlockMatrix(const BlockMatrix& other);

    void resize(unsigned int numRows, unsigned int numCols);

    unsigned int getNumberRows();

    unsigned int getNumberCols();

    MatrixEpetraPtr& block(unsigned int row, unsigned int col);

    MatrixEpetraPtr block(unsigned int row, unsigned int col) const;

    void add(const BlockMatrix& other);

    // void multiply(const BlockMatrix& other, BlockMatrix& result);

    VectorEpetra operator*(const VectorEpetra& vector);

    void copyBlock(unsigned int rows, unsigned int cols,
                   MatrixEpetraPtr matrix);

    MapPtr rangeMap(unsigned int row, unsigned int col) const;

    MapPtr domainMap(unsigned int row, unsigned int col) const;

    BlockMatrix& operator*=(const double& coeff);

    Grid getGrid();

    void spy(std::string prefix = "");

    void singleNorms1();

    void setMaps(std::vector<MapPtr> rangeMaps, std::vector<MapPtr> domainMaps);

    void printPattern();

private:
    Grid                 M_grid;
    GridEpetra           M_gridEpetra;
    unsigned int         M_rows;
    unsigned int         M_cols;
    std::vector<MapPtr>  M_rangeMaps;
    std::vector<MapPtr>  M_domainMaps;
    MapPtr               M_globalDomainMap;
    MapPtr               M_globalRangeMap;
};

}  // namespace RedMA

#endif  // BLOCKMATRIX_HPP
