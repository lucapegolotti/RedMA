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

#include <rb/reduced_basis/rbSolver/ProperOrthogonalDecomposition.hpp>

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

    std::vector<SHP(VECTOREPETRA)> getBasis(const unsigned int& index);

    std::vector<SHP(VECTOREPETRA)> getFullBasis(const unsigned int& index) {return M_bases[index];}

    std::vector<SHP(VECTOREPETRA)> getEnrichedBasis(const unsigned int& index);

    std::vector<SHP(VECTOREPETRA)> getPrimalSupremizers(const unsigned int& i, const unsigned int& j) {return M_primalSupremizers(i,j);}

    DenseMatrix matrixProject(MatrixEp matrix, unsigned int basisIndexRow,
                                               unsigned int basisIndexCol);

    BlockMatrix<DenseMatrix> leftProject(BlockMatrix<MatrixEp> matrix);

    DenseMatrix leftProject(MatrixEp matrix, unsigned int basisIndex);

    BlockVector<DenseVector> leftProject(BlockVector<VectorEp> vector);

    DenseVector leftProject(VectorEp vector, unsigned int basisIndex);

    BlockVector<DenseVector> projectOnLagrangeSpace(BlockVector<VectorEp> vector);

    BlockMatrix<DenseMatrix> rightProject(BlockMatrix<MatrixEp> matrix);

    DenseMatrix rightProject(MatrixEp matrix, unsigned int basisIndex);

    SHP(VECTOREPETRA) reconstructFEFunction(DenseVector rbSolution, unsigned int index);

    void addPrimalSupremizer(SHP(VECTOREPETRA) supremizer,
                             const unsigned int& fieldToAugment,
                             const unsigned int& fieldConstraint);

    void addDualSupremizer(SHP(VECTOREPETRA) supremizer,
                           const unsigned int& fieldToAugment);

    inline unsigned int getSizeBasis(const unsigned int& index) {return getBasis(index).size();}

    inline unsigned int getFullSizeBasis(const unsigned int& index) {return M_bases[index].size();}

    inline unsigned int getSizeEnrichedBasis(const unsigned int& index) {return getEnrichedBasis(index).size();}

    void normalizeBasis(const unsigned int& index, SHP(MATRIXEPETRA) normMatrix);

    void computeOnlineNumberBasisFunctions(unsigned int index);

    // indices of the vectors from the ith enriched basis to keep
    std::vector<unsigned int> getSelectors(unsigned int index);

    void print();

    void setBasisFunctions(std::vector<SHP(VECTOREPETRA)> basisFunctions,
                           unsigned int index) {M_bases[index] = basisFunctions;}

    inline unsigned int getNumFields() {return M_numFields;}

    inline bool hasSupremizers(const unsigned int& index)
    {
        for (unsigned int j = 0; j < M_numFields; j++)
        {
            if (M_primalSupremizers(index,j).size() > 0)
                return true;
        }

        if (M_dualSupremizers[index].size() > 0)
            return true;

        return false;
    }

private:
    void addVectorsFromFile(std::string filename,
                            std::vector<SHP(VECTOREPETRA)>& vectors,
                            const unsigned int& indexField);

    unsigned int                                    M_numFields;
    DataContainer                                   M_data;
    EPETRACOMM                                      M_comm;
    std::vector<std::vector<SHP(VECTOREPETRA)>>     M_bases;
    std::vector<MatrixEp>                           M_enrichedBasesMatrices;
    std::vector<MatrixEp>                           M_enrichedBasesMatricesTransposed;
    // this is a grid because the row indicates the field to be augmented (velocity)
    // and the column indicates the constraining field (pressure)
    GridVectors                                     M_primalSupremizers;
    std::vector<std::vector<SHP(VECTOREPETRA)>>     M_dualSupremizers;
    std::string                                     M_meshName;
    std::string                                     M_path;
    std::vector<std::vector<double>>                M_svs;
    std::vector<SHP(FESPACE)>                       M_fespaces;
    // max of tolerance on the individual fields
    double                                          M_onlineTol;
    // online tolerance on the intividual fields
    std::vector<double>                             M_onlineTols;
    std::vector<unsigned int>                       M_NsOnline;
};

}  // namespace RedMA

#endif  // RBBASES_HPP
