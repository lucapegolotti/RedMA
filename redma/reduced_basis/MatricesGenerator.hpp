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

#ifndef MATRICESGENERATOR_HPP
#define MATRICESGENERATOR_HPP

#include <redma/RedMA.hpp>
#include <redma/utils/PrintLog.hpp>
#include <redma/problem/DataContainer.hpp>
#include <redma/problem/GlobalProblem.hpp>
#include <redma/assemblers/AssemblerFactory.hpp>
#include <redma/assemblers/coupling/InterfaceAssembler.hpp>

#include <redma/geometry/TreeStructure.hpp>
#include <redma/geometry/building_blocks/BuildingBlock.hpp>
#include <redma/geometry/GeometryPrinter.hpp>
#include <redma/geometry/building_blocks/Tube.hpp>
#include <redma/geometry/building_blocks/BifurcationSymmetric.hpp>

#include <redma/reduced_basis/RBBases.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

#include <redma/RedMA.hpp>

#include <fstream>

namespace RedMA
{

/*! \brief Class for generating the finite element matrices necessary to the
 * offline phase of the reduced basis method.
 */
class MatricesGenerator
{
    typedef aAssembler                                      AssemblerType;
    typedef std::vector<std::vector<shp<VECTOREPETRA>>>     VectorFunctions;
    typedef std::pair<shp<AssemblerType>, VectorFunctions>  AssemblerSnapshotPair;
public:
    /*! \brief Constructor.
     *
     * \param data A DataContainer.
     * \param comm The MPI Communicator.
     */
    MatricesGenerator(const DataContainer& data,
                      EPETRACOMM comm);

    /*! \brief Generate matrices related to problem (the one specified in the
     * geometry file).
     */
    void generate();

private:
    void createDefaultAssemblers();

    shp<TreeNode> generateDefaultTreeNode(const std::string& nameMesh);

    shp<TreeNode> generateDefaultTube(const std::string& nameMesh);

    shp<TreeNode> generateDefaultSymmetricBifurcation(const std::string& nameMesh);

    DataContainer                                       M_data;
    EPETRACOMM                                          M_comm;
    std::map<std::string, AssemblerSnapshotPair>        M_meshASPairMap;
    std::map<std::string, shp<RBBases>>                 M_bases;
};

}  // namespace RedMA

#endif  // MATRICESGENERATOR_HPP
