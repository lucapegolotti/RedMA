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

#ifndef BASISGENERATOR_HPP
#define BASISGENERATOR_HPP

#include <redma/RedMA.hpp>
#include <redma/utils/PrintLog.hpp>
#include <redma/solver/problem/DataContainer.hpp>
#include <redma/solver/problem/ProblemFEM.hpp>
#include <redma/solver/assemblers/AssemblerFactory.hpp>
#include <redma/solver/assemblers/InterfaceAssembler.hpp>

#include <redma/geometry/TreeStructure.hpp>
#include <redma/geometry/BuildingBlock.hpp>
#include <redma/geometry/GeometryPrinter.hpp>
#include <redma/geometry/Tube.hpp>
#include <redma/geometry/BifurcationSymmetric.hpp>

#include <rb/reduced_basis/rbSolver/ProperOrthogonalDecomposition.hpp>
#include <redma/reduced_basis/RBBases.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

#include <boost/filesystem.hpp>

#include <fstream>

namespace RedMA
{

class BasisGenerator
{
    typedef aAssembler<FEVECTOR COMMA FEMATRIX>             AssemblerType;
    typedef std::vector<std::vector<SHP(VECTOREPETRA)>>     VectorFunctions;
    typedef std::pair<SHP(AssemblerType), VectorFunctions>  AssemblerSnapshotPair;
public:
    BasisGenerator(const DataContainer& data, EPETRACOMM comm);

    void generateBasis();

    // use this method when you want to only dump the matrices in order to compute
    // for example the basis with matlab
    void generateMatricesOnly();

private:
    void createDefaultAssemblers();

    void parseFiles();

    void performPOD();

    void dumpBasis();

    void addSupremizers();

    void orthonormalize();

    SHP(TreeNode) generateDefaultTreeNode(const std::string& nameMesh);

    SHP(TreeNode) generateDefaultTube(const std::string& nameMesh);

    SHP(TreeNode) generateDefaultSymmetricBifurcation(const std::string& nameMesh);

    void parseParameterSnapshots(const std::string& paramDir);

    void addSnapshotsFromFile(const std::string& snapshotsFile,
                              std::vector<SHP(VECTOREPETRA)>& snapshots,
                              SHP(FESPACE) fespace);

    LifeV::LinearSolver setupLinearSolver(MatrixEp matrix);

    DataContainer                                       M_data;
    EPETRACOMM                                          M_comm;
    std::map<std::string, AssemblerSnapshotPair>        M_meshASPairMap;
    std::map<std::string, SHP(RBBases)>                 M_bases;
};

}  // namespace RedMA

#endif  // BASISGENERATOR_HPP
