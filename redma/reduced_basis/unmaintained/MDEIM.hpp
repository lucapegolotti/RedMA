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
#include <redma/problem/DataContainer.hpp>

#include <redma/array/SparseMatrix.hpp>

//  #include <rb/reduced_basis/rbSolver/ProperOrthogonalDecomposition.hpp>
//  #include <rb/reduced_basis/util/EpetraArrayUtils.hpp>
#include <redma/reduced_basis/MDEIMStructure.hpp>

#include <redma/RedMA.hpp>

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

    void prepareOnline();

    void checkOnline(FEMATRIX reducedMatrix, FEMATRIX fullMatrix);

    void initialize(FEMATRIX matrix);

    void setFESpace(shp<FESPACE> fespace);

    void checkOnSnapshots();

    shp<MDEIMStructure>& getMDEIMStructure() {return M_structure;}

    void dumpMDEIM(std::string dir);

    void projectMDEIM(std::vector<shp<VECTOREPETRA>> leftBasis,
                      std::vector<shp<VECTOREPETRA>> rightBasis);

    void setDomainMap(shp<MAPEPETRA> domainMap) {M_domainMap = domainMap;}

    void setRangeMap(shp<MAPEPETRA> rangeMap) {M_rangeMap = rangeMap;}

    void loadMDEIM(std::string pathdir);

    FEMATRIX assembleMatrix(FEMATRIX reducedMatrix);

    RBMATRIX assembleProjectedMatrix(FEMATRIX reducedMatrix);

private:

    shp<VECTOREPETRA> vectorizeMatrix(FEMATRIX matrix);

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

    void loadProjectedBasis(std::string filename);

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

    FEMATRIX                                    M_firstSnapshot;
    std::vector<shp<VECTOREPETRA>>              M_snapshotsVectorized;
    std::vector<shp<VECTOREPETRA>>              M_basis;
    std::vector<shp<DENSEVECTOR>>               M_basisProjected;
    shp<MDEIMStructure>                         M_structure;
    EPETRACOMM                                  M_comm;
    DataContainer                               M_data;
    shp<FESPACE>                                M_fespace;
    bool                                        M_isInitialized;
    shp<MAPEPETRA>                              M_domainMap;
    shp<MAPEPETRA>                              M_rangeMap;
};

}  // namespace RedMA

#endif  // MDEIM_HPP
