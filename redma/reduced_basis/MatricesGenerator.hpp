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

#include "aMatricesGenerator.hpp"

namespace RedMA
{

/*! \brief Class for generating the finite element matrices necessary to the
 * offline phase of the reduced basis method.
 */
class MatricesGenerator : public aMatricesGenerator
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
    virtual void generate() override;

protected:
    virtual void createAssemblers() override;

    shp<TreeNode> generateDefaultTreeNode(const std::string& nameMesh);

    shp<TreeNode> generateDefaultTube(const std::string& nameMesh);

    shp<TreeNode> generateDefaultSymmetricBifurcation(const std::string& nameMesh);

    shp<TreeNode> generateDefaultAortaBifurcation0(const std::string& nameMesh);

    shp<TreeNode> generateDefaultAortaBifurcation1(const std::string& nameMesh);

    shp<TreeNode> generateDefaultBypass(const std::string& nameMesh);

    std::map<std::string, AssemblerSnapshotPair>        M_meshASPairMap;
    std::map<std::string, shp<RBBases>>                 M_bases;
};

}  // namespace RedMA

#endif  // MATRICESGENERATOR_HPP
