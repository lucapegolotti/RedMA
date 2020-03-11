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

#ifndef RBBASES_HPP
#define RBBASES_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/problem/DataContainer.hpp>

#include <redma/solver/array/BlockMatrix.hpp>

#include <boost/filesystem.hpp>

namespace RedMA
{

class RBBases
{
    typedef boost::numeric::ublas::matrix<std::vector<SHP(VECTOREPETRA)>> GridVectors;

public:
    RBBases(const DataContainer& data, EPETRACOMM comm);

    void setNumberOfFields(const unsigned int& nfields);

    void setPath(std::string path);

    void loadSingularValues();

    void loadBases();

    void dump();

    void setFESpace(SHP(FESPACE) fespace, const unsigned int& indexbasis);

    std::vector<SHP(VECTOREPETRA)>& getBasis(const unsigned int& index, double tol = 0);

    std::vector<SHP(VECTOREPETRA)> getEnrichedBasis(const unsigned int& index, double tol = 0);

    std::vector<SHP(VECTOREPETRA)> getPrimalSupremizers(const unsigned int& i, const unsigned int& j) {return M_primalSupremizers(i,j);}

    BlockMatrix<DenseMatrix> leftProject(BlockMatrix<MatrixEp> matrix);

    DenseMatrix leftProject(MatrixEp matrix, unsigned int basisIndex);

    BlockVector<DenseVector> leftProject(BlockVector<VectorEp> vector);

    DenseVector leftProject(VectorEp vector, unsigned int basisIndex);

    BlockMatrix<DenseMatrix> rightProject(BlockMatrix<MatrixEp> matrix);

    DenseMatrix rightProject(MatrixEp matrix, unsigned int basisIndex);

    void addPrimalSupremizer(SHP(VECTOREPETRA) supremizer,
                             const unsigned int& fieldToAugment,
                             const unsigned int& fieldConstraint);

    void addDualSupremizer(SHP(VECTOREPETRA) supremizer,
                           const unsigned int& fieldToAugment);

    inline unsigned int getSizeBasis(const unsigned int& index) {return M_bases[index].size();}

    inline unsigned int getSizeEnrichedBasis(const unsigned int& index) {return getEnrichedBasis(index).size();}

private:
    void addVectorsFromFile(std::string filename,
                            std::vector<SHP(VECTOREPETRA)>& vectors,
                            const unsigned int& indexField);

    unsigned int                                    M_numFields;
    DataContainer                                   M_data;
    EPETRACOMM                                      M_comm;
    std::vector<std::vector<SHP(VECTOREPETRA)>>     M_bases;
    // this is a grid because the row indicates the field to be augmented (velocity)
    // and the column indicates the constraining field (pressure)
    GridVectors                                     M_primalSupremizers;
    std::vector<std::vector<SHP(VECTOREPETRA)>>     M_dualSupremizers;
    std::string                                     M_meshName;
    std::string                                     M_path;
    std::vector<std::vector<double>>                M_svs;
    std::vector<SHP(FESPACE)>                       M_fespaces;
};

}  // namespace RedMA

#endif  // RBBASES_HPP
