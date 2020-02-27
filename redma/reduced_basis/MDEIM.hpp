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
#include <redma/reduced_basis/MDEIMStructure.hpp>

namespace RedMA
{

class MDEIM
{
public:
    MDEIM();

    void setDataContainer(const DataContainer& dataContainer);

    void setComm(EPETRACOMM comm);

    void addSnapshot(MatrixEp newSnapshot);

    void performMDEIM();

    void prepareOnline(MatrixEp matrix);

    void checkOnline(MatrixEp reducedMatrix, MatrixEp fullMatrix);

    void initialize(MatrixEp matrix);

    void setFESpace(SHP(FESPACE) fespace);

    void checkOnSnapshots();

    SHP(MDEIMStructure)& getMDEIMStructure() {return M_structure;}

    void dumpMDEIM(std::string dir);

private:

    void vectorizeSnapshots();

    SHP(VECTOREPETRA) vectorizeMatrix(MatrixEp matrix);

    void performPOD();

    void pickMagicPoints();

    void computeInterpolationVectorOffline(VECTOREPETRA& vector,
                                           Epetra_SerialDenseVector& interpolationCoefficients,
                                           Epetra_SerialDenseSolver& solver);

    void computeFeInterpolation(Epetra_SerialDenseVector& interpolationCoefficients,
                                VECTOREPETRA& vector);

    void buildReducedMesh();

    void identifyReducedNodes();

    void identifyReducedElements();

    void computeInterpolationRhsOnline(Epetra_SerialDenseVector& interpVector,
                                       MatrixEp reducedMat);

    void computeInterpolationVectorOnline(Epetra_SerialDenseVector& interpVector,
                                          MatrixEp reducedMat);

    void reconstructMatrixFromVectorizedForm(VECTOREPETRA& vectorizedAh,
                                             MATRIXEPETRA& Ah);

    std::vector<MatrixEp>                       M_snapshots;
    std::vector<SHP(VECTOREPETRA)>              M_snapshotsVectorized;
    std::vector<SHP(VECTOREPETRA)>              M_basis;
    SHP(MDEIMStructure)                         M_structure;
    EPETRACOMM                                  M_comm;
    DataContainer                               M_data;
    SHP(FESPACE)                                M_fespace;
    bool                                        M_isInitialized;
};

}  // namespace RedMA

#endif  // MDEIM_HPP
