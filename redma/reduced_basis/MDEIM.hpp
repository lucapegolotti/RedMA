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
#include <redma/solver/assemblers/AssemblerFactory.hpp>

#include <rb/reduced_basis/rbSolver/ProperOrthogonalDecomposition.hpp>
#include <rb/reduced_basis/util/EpetraArrayUtils.hpp>

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
            delete[] myRowMatrixEntriesOfMagicPoints;
            delete[] myColMatrixEntriesOfMagicPoints;
            delete[] rowLocalReducedIndeces;
            delete[] globalReducedNodes;
            delete[] myGlobalReducedNodes;
            delete[] reducedElements;

            for (unsigned int i = 0; i < numMyRows; i++)
                delete[] columnIndeces[i];

            delete[] columnIndeces;
        }
    }

    int                             N;
    bool                            allocated;
    int                             numGlobalNonzeros;
    int                             numMyNonzeros;
    int                             numMyRows;
    int*                            numMyEntries;
    int**                           columnIndeces;
    int*                            partialSumMyEntries;
    int                             numMyLocalMagicPoints;
    std::vector<int>                myLocalMagicPoints;
    std::vector<int>                localIndecesMagicPoints;
    std::vector<int>                magicPointsProcOwner;
    std::vector<int>                globalIndecesMagicPoints;
    int*                            myRowMatrixEntriesOfMagicPoints;
    int*                            myColMatrixEntriesOfMagicPoints;
    int*                            rowLocalReducedIndeces;
    int*                            globalReducedNodes;
    int*                            myGlobalReducedNodes;
    int                             numGlobalReducedNodes;
    int                             numMyGlobalReducedNodes;
    int                             numReducedElements;
    unsigned int*                   reducedElements;
    Epetra_SerialDenseMatrix        Qj;
    SHP(MAPEPETRA)                  vectorMap;
};

class MDEIM
{
    typedef boost::numeric::ublas::matrix<SHP(SingleMDEIMStructure)>      GridStructures;
    typedef boost::numeric::ublas::matrix<std::vector<SHP(VECTOREPETRA)>> GridVectors;

public:
    MDEIM();

    void setDataContainer(const DataContainer& dataContainer);

    void setComm(EPETRACOMM comm);

    void setAssembler(SHP(aAssembler<FEVECTOR COMMA FEMATRIX>));

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

    void pickMagicPoints(const unsigned int& i, const unsigned int& j);

    void computeInterpolationVectorOffline(VECTOREPETRA& vector,
                                           Epetra_SerialDenseVector& interpolationCoefficients,
                                           Epetra_SerialDenseSolver& solver,
                                           SHP(SingleMDEIMStructure) mstruct);

    void computeFeInterpolation(const unsigned int& i, const unsigned int& j,
                                Epetra_SerialDenseVector& interpolationCoefficients,
                                VECTOREPETRA& vector);

    void buildReducedMesh(const unsigned int& i, const unsigned int& j);

    void identifyReducedNodes(const unsigned int& i, const unsigned int& j);

    void identifyReducedElements(const unsigned int& i, const unsigned int& j);

    std::vector<BlockMatrix<MatrixEp>>          M_snapshots;
    GridVectors                                 M_snapshotsVectorized;
    GridVectors                                 M_bases;
    GridStructures                              M_structures;
    EPETRACOMM                                  M_comm;
    DataContainer                               M_data;
    SHP(aAssembler<FEVECTOR COMMA FEMATRIX>)    M_assembler;
    // rows and cols of the block structures
    unsigned int                                M_nRows;
    unsigned int                                M_nCols;
};

}  // namespace RedMA

#endif  // MDEIM_HPP
