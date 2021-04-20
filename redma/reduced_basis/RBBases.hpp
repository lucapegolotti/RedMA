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
#include <redma/problem/DataContainer.hpp>

#include <redma/array/BlockMatrix.hpp>

// #include <rb/reduced_basis/rbSolver/ProperOrthogonalDecomposition.hpp>

#include <redma/RedMA.hpp>

namespace RedMA
{

class RBBases
{
    typedef boost::numeric::ublas::matrix<std::vector<shp<VECTOREPETRA>>> GridVectors;

public:
    RBBases(const DataContainer& data, EPETRACOMM comm);

    void setNumberOfFields(const unsigned int& nfields);

    void setPath(std::string path);

    void loadSingularValues();

    void loadBases();

    void dump();

    void setFESpace(shp<FESPACE> fespace, const unsigned int& indexbasis);

    std::vector<shp<VECTOREPETRA>> getBasis(const unsigned int& index);

    std::vector<shp<VECTOREPETRA>> getFullBasis(const unsigned int& index) {return M_bases[index];}

    std::vector<shp<VECTOREPETRA>> getEnrichedBasis(const unsigned int& index, unsigned int ID);

    std::vector<shp<VECTOREPETRA>> getPrimalSupremizers(const unsigned int& i, const unsigned int& j) {return M_primalSupremizers(i,j);}

    shp<aMatrix> matrixProject(shp<aMatrix> matrix, unsigned int basisIndexRow,
                               unsigned int basisIndexCol, unsigned int ID);

    shp<BlockMatrix> leftProject(shp<BlockMatrix> matrix, unsigned int ID);

    shp<DenseMatrix> leftProject(shp<SparseMatrix> matrix, unsigned int basisIndex, unsigned int ID);

    shp<BlockVector> leftProject(shp<BlockVector> vector, unsigned int ID);

    shp<DenseVector> leftProject(shp<DistributedVector> vector, unsigned int basisIndex, unsigned int ID);

    shp<BlockVector> projectOnLagrangeSpace(shp<BlockVector> vector);

    shp<BlockMatrix> rightProject(shp<BlockMatrix> matrix, unsigned int ID);

    shp<DenseMatrix> rightProject(shp<SparseMatrix> matrix, unsigned int basisIndex, unsigned int ID);

    shp<VECTOREPETRA> reconstructFEFunction(shp<aVector> rbSolution, unsigned int index,
                                            unsigned int ID);

    void addPrimalSupremizer(shp<VECTOREPETRA> supremizer,
                             const unsigned int& fieldToAugment,
                             const unsigned int& fieldConstraint);

    void addDualSupremizer(shp<VECTOREPETRA> supremizer,
                           const unsigned int& fieldToAugment);

    inline unsigned int getSizeBasis(const unsigned int& index) {return getBasis(index).size();}

    inline unsigned int getFullSizeBasis(const unsigned int& index) {return M_bases[index].size();}

    inline unsigned int getSizeEnrichedBasis(const unsigned int& index) {return getEnrichedBasis(index,0).size();}

    void normalizeBasis(const unsigned int& index, shp<MATRIXEPETRA> normMatrix);

    void computeOnlineNumberBasisFunctions(unsigned int index);

    // indices of the vectors from the ith enriched basis to keep
    std::vector<unsigned int> getSelectors(unsigned int index);

    void print();

    void scaleBasisWithPiola(unsigned int index, unsigned int ID,
                             std::function<void(shp<VECTOREPETRA>)> transform);

    void setBasisFunctions(std::vector<shp<VECTOREPETRA>> basisFunctions,
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

    shp<SparseMatrix> getEnrichedBasisMatrices(const unsigned int& index, const unsigned int& ID,
                                               bool transpose);

private:
    void addVectorsFromFile(std::string filename,
                            std::vector<shp<VECTOREPETRA>>& vectors,
                            const unsigned int& indexField,
                            int Nmax = -1);

    unsigned int                                                      M_numFields;
    DataContainer                                                     M_data;
    EPETRACOMM                                                        M_comm;
    std::vector<std::vector<shp<VECTOREPETRA>>>                       M_bases;
    std::vector<shp<SparseMatrix>>                                    M_enrichedBasesMatrices;
    std::vector<shp<SparseMatrix>>                                    M_enrichedBasesMatricesTransposed;
    // this is a grid because the row indicates the field to be augmented (velocity)
    // and the column indicates the constraining field (pressure)
    GridVectors                                                       M_primalSupremizers;
    std::vector<std::vector<shp<VECTOREPETRA>>>                       M_dualSupremizers;
    std::string                                                       M_meshName;
    std::string                                                       M_path;
    std::vector<std::vector<double>>                                  M_svs;
    std::vector<shp<FESPACE>>                                         M_fespaces;
    bool                                                              M_fespacesAreSet;
    // max of tolerance on the individual fields
    double                                                            M_onlineTol;
    // online tolerance on the individual fields
    std::vector<double>                                               M_onlineTols;
    std::vector<unsigned int>                                         M_NsOnline;
    // these are needed only when we use a different basis in every block (e.g.
    // when we rescale with piola)
    std::map<unsigned int,std::map<unsigned int,
                          std::vector<shp<VECTOREPETRA>>>>            M_enrichedBasesMap;
    std::map<unsigned int,std::map<unsigned int,shp<SparseMatrix>>>   M_enrichedBasesMatricesMap;
    std::map<unsigned int,std::map<unsigned int,shp<SparseMatrix>>>   M_enrichedBasesMatricesTransposedMap;
};

}  // namespace RedMA

#endif  // RBBASES_HPP
