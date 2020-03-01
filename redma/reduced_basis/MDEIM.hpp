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

    void addSnapshot(FEMATRIX newSnapshot);

    void performMDEIM(std::string outdir);

    void prepareOnline(FEMATRIX matrix);

    void checkOnline(FEMATRIX reducedMatrix, FEMATRIX fullMatrix);

    void initialize(FEMATRIX matrix);

    void setFESpace(SHP(FESPACE) fespace);

    void checkOnSnapshots();

    SHP(MDEIMStructure)& getMDEIMStructure() {return M_structure;}

    void dumpMDEIM(std::string dir);

    void projectMDEIM(std::vector<SHP(VECTOREPETRA)> leftBasis,
                      std::vector<SHP(VECTOREPETRA)> rightBasis);

    void setDomainMap(SHP(MAPEPETRA) domainMap) {M_domainMap = domainMap;}

    void setRangeMap(SHP(MAPEPETRA) rangeMap) {M_rangeMap = rangeMap;}

    void loadMDEIM(std::string pathdir);

    FEMATRIX assembleMatrix(FEMATRIX reducedMatrix);

    RBMATRIX assembleProjectedMatrix(FEMATRIX reducedMatrix);

private:

    void vectorizeSnapshots();

    SHP(VECTOREPETRA) vectorizeMatrix(FEMATRIX matrix);

    void performPOD(std::string outdir);

    void pickMagicPoints();

    void computeInterpolationVectorOffline(VECTOREPETRA& vector,
                                           DENSEVECTOR& interpolationCoefficients,
                                           DENSESOLVER& solver);

    void computeFeInterpolation(DENSEVECTOR& interpolationCoefficients,
                                VECTOREPETRA& vector);

    void computeProjectedInterpolation(DENSEVECTOR& interpolationCoefficients,
                                       DENSEVECTOR& vector);

    void loadBasis(std::string filename);

    void buildReducedMesh();

    void identifyReducedNodes();

    void identifyReducedElements();

    void computeInterpolationRhsOnline(DENSEVECTOR& interpVector,
                                       FEMATRIX reducedMat);

    void computeInterpolationVectorOnline(DENSEVECTOR& interpVector,
                                          FEMATRIX reducedMat);

    void reconstructMatrixFromVectorizedForm(VECTOREPETRA& vectorizedAh,
                                             MATRIXEPETRA& Ah);

    void reconstructMatrixFromVectorizedForm(DENSEVECTOR& vectorizedAn,
                                             DENSEMATRIX& An);

    void dumpBasis(std::string dir);

    void dumpProjectedBasis(std::string dir);

    std::vector<FEMATRIX>                       M_snapshots;
    std::vector<SHP(VECTOREPETRA)>              M_snapshotsVectorized;
    std::vector<SHP(VECTOREPETRA)>              M_basis;
    std::vector<SHP(DENSEVECTOR)>               M_basisProjected;
    SHP(MDEIMStructure)                         M_structure;
    EPETRACOMM                                  M_comm;
    DataContainer                               M_data;
    SHP(FESPACE)                                M_fespace;
    bool                                        M_isInitialized;
    SHP(MAPEPETRA)                              M_domainMap;
    SHP(MAPEPETRA)                              M_rangeMap;
};

}  // namespace RedMA

#endif  // MDEIM_HPP
