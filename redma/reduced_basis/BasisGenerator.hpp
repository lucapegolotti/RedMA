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

#include <redma/geometry/TreeStructure.hpp>
#include <redma/geometry/GeometryPrinter.hpp>
#include <redma/geometry/Tube.hpp>

#include <rb/reduced_basis/rbSolver/ProperOrthogonalDecomposition.hpp>

#include <boost/filesystem.hpp>

#include <fstream>

namespace RedMA
{

class BasisGenerator
{
    typedef aAssembler<FEVECTOR COMMA FEMATRIX>           AssemblerType;
    typedef std::vector<std::vector<SHP(VECTOREPETRA)>>   Snapshots;
    typedef std::pair<SHP(AssemblerType), Snapshots>      AssemblerSnapshotPair;
public:
    BasisGenerator(const DataContainer& data, EPETRACOMM comm);

    void generateBasis();


private:
    void parseFiles();

    void performPOD();

    SHP(TreeNode) generateDefaultTreeNode(const std::string& nameMesh);

    SHP(TreeNode) generateDefaultTube(const std::string& nameMesh);

    void parseParameterSnapshots(const std::string& paramDir);

    void addSnapshotsFromFile(const std::string& snapshotsFile,
                              std::vector<SHP(VECTOREPETRA)>& snapshots,
                              SHP(FESPACE) fespace);

    DataContainer                                       M_data;
    EPETRACOMM                                          M_comm;
    std::map<std::string, AssemblerSnapshotPair>        M_meshASPairMap;
};

}  // namespace RedMA

#endif  // BASISGENERATOR_HPP
