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

#ifndef MDEIM_HPP
#define MDEIM_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/problem/DataContainer.hpp>
#include <redma/solver/problem/ProblemFEM.hpp>

#include <rb/reduced_basis/rbSolver/ProperOrthogonalDecomposition.hpp>

namespace RedMA
{

struct SingleMDEIMStructure
{
    SingleMDEIMStructure() :
      allocated(false)
    {

    }

    ~SingleMDEIMStructure()
    {
        if (allocated)
        {
            delete[] numMyEntries;
            delete[] partialSumMyEntries;

            for (unsigned int i = 0; i < numMyRows; i++)
                delete[] columnIndeces[i];

            delete[] columnIndeces;
        }
    }

    bool           allocated;
    int            numGlobalNonzeros;
    int            numMyNonzeros;
    int            numMyRows;
    int*           numMyEntries;
    int**          columnIndeces;
    int*           partialSumMyEntries;
    SHP(MAPEPETRA) vectorMap;
};

class MDEIM
{
    typedef boost::numeric::ublas::matrix<SHP(SingleMDEIMStructure)>      GridStructures;
    typedef boost::numeric::ublas::matrix<std::vector<SHP(VECTOREPETRA)>> GridVectors;

public:
    MDEIM();

    void setDataContainer(const DataContainer& dataContainer);

    void setComm(EPETRACOMM comm);

    void addSnapshot(BlockMatrix<MatrixEp> newSnapshot);

    void performMDEIM();

private:

    void initializeMDEIMStructures(BlockMatrix<MatrixEp> matrix);

    void initializeSingleMDEIMStructure(const unsigned int& i,
                                        const unsigned int& j,
                                        MatrixEp matrix);

    void vectorizeSnapshots();

    SHP(VECTOREPETRA) vectorizeMatrix(const unsigned int& i,
                                      const unsigned int& j,
                                      SHP(MATRIXEPETRA) matrix);

    void performPOD();

    std::vector<BlockMatrix<MatrixEp>>          M_snapshots;
    GridVectors                                 M_snapshotsVectorized;
    GridVectors                                 M_bases;
    GridStructures                              M_structures;
    EPETRACOMM                                  M_comm;
    DataContainer                               M_data;
    // rows and cols of the block structures
    unsigned int                                M_nRows;
    unsigned int                                M_nCols;
};

}  // namespace RedMA

#endif  // MDEIM_HPP
