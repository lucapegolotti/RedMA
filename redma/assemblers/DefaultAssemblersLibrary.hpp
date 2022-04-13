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

#ifndef DEFAULTASSEMBLERSLIBRARY_HPP
#define DEFAULTASSEMBLERSLIBRARY_HPP

#include <redma/RedMA.hpp>
#include <redma/problem/DataContainer.hpp>

#include <redma/geometry/TreeStructure.hpp>

namespace RedMA
{

class aAssembler;

shp<aAssembler>
AssemblerFactory(const DataContainer& data, shp<TreeNode> treeNode);

/*! \brief Class to handle a set of default (i.e., defined on reference
 *   building blocks) assemblers.
 */
class DefaultAssemblersLibrary
{
    typedef aAssembler   AssemblerType;
public:
    /*! \brief Constructor.
     *
     * \param data A DataContainer.
     * \param meshes Set of mesh names.
     * \param comm MPI Communicator.
     */
    DefaultAssemblersLibrary(const DataContainer& data,
                             const std::set<std::string>& meshes,
                             EPETRACOMM comm);

    /*! \brief Generate a tree node for a given mesh.
     *
     * \param nameMesh Name of the mesh.
     * \return Shared pointer to the TreeNode.
     */
    shp<TreeNode> generateDefaultTreeNode(const std::string& nameMesh);

    /*! \brief Generate a tree node corresponding to a tube.
     *
     * \param nameMesh Name of the mesh.
     * \return Shared pointer to the TreeNode.
     */
    shp<TreeNode> generateDefaultTube(const std::string& nameMesh);

    /*! \brief Generate a tree node corresponding to a symmetric bifurcation.
     *
     * \param nameMesh Name of the mesh.
     * \return Shared pointer to the TreeNode.
     */
    shp<TreeNode> generateDefaultSymmetricBifurcation(const std::string& nameMesh);

    /*! \brief Generate a tree node corresponding to a femoropopliteal bypass.
     *
     * \param nameMesh Name of the mesh.
     * \return Shared pointer to the TreeNode.
     */
    shp<TreeNode> generateDefaultBypass(const std::string& nameMesh);

    /*! \brief Getter for a given default assembler.
     *
     * \param nameMesh Name of the mesh.
     * \return Shared pointer desired assembler.
     */
    shp<AssemblerType> getDefaultAssembler(const std::string& namemesh);

private:
    DataContainer                                   M_data;
    std::map<std::string, shp<AssemblerType>>       M_assemblersMap;
    EPETRACOMM                                      M_comm;
    unsigned int                                    M_count;
};

}

#endif // DEFAULTASSEMBLERSLIBRARY_HPP
