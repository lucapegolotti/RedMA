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
#include <redma/reduced_basis/SingleMDEIMStructure.hpp>

namespace RedMA
{

class MDEIM
{
    typedef boost::numeric::ublas::matrix<std::vector<SHP(VECTOREPETRA)>> GridVectors;

public:
    MDEIM();

    void setDataContainer(const DataContainer& dataContainer);

    void setComm(EPETRACOMM comm);

    void setAssembler(SHP(aAssembler<FEVECTOR COMMA FEMATRIX>));

    void addSnapshot(BlockMatrix<MatrixEp> newSnapshot);

    void performMDEIM();

    void prepareOnline();

    void checkOnline();

    void checkOnSnapshots();

    inline void setMatrixIndex(const unsigned int& index) {M_matIndex = index;}

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

    void computeInterpolationRhsOnline(const unsigned int& i,
                                       const unsigned int& j,
                                       Epetra_SerialDenseVector& interpVector);

    void prepareOnline(const unsigned int& i, const unsigned int& j, BlockMatrix<FEMATRIX> mat);

    void checkOnline(const unsigned int& i, const unsigned int& j);

    void computeInterpolationVectorOnline(const unsigned int& i,
                                          const unsigned int& j,
                                          Epetra_SerialDenseVector& interpVector);

    void reconstructMatrixFromVectorizedForm(const unsigned int& i,
                                             const unsigned int& j,
                                             VECTOREPETRA& vectorizedAh,
                                             MATRIXEPETRA& Ah);

    std::vector<BlockMatrix<MatrixEp>>          M_snapshots;
    GridVectors                                 M_snapshotsVectorized;
    GridVectors                                 M_bases;
    GridMDEIMStructures                         M_structures;
    EPETRACOMM                                  M_comm;
    DataContainer                               M_data;
    SHP(aAssembler<FEVECTOR COMMA FEMATRIX>)    M_assembler;
    unsigned int                                M_matIndex;
    // rows and cols of the block structures
    unsigned int                                M_nRows;
    unsigned int                                M_nCols;
};

}  // namespace RedMA

#endif  // MDEIM_HPP
