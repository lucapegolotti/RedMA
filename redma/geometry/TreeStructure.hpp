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

#ifndef TREESTRUCTURE_HPP
#define TREESTRUCTURE_HPP

#include <redma/geometry/building_blocks/BuildingBlock.hpp>
#include <redma/geometry/building_blocks/Tube.hpp>
#include <redma/geometry/building_blocks/BifurcationSymmetric.hpp>

#include <Epetra_SerialComm.h>
#include <Epetra_MpiComm.h>

#include <memory>
#include <map>
#include <set>
#include <queue>

#include <ctime>
#include <cstdlib>

namespace RedMA
{

/// A node containing a geometric building block.
class TreeNode
{
private:
    typedef LifeV::VectorSmall<3>          Vector3D;
public:

    /*! \brief Constructor.
     *
     * \param block Shared pointer to the building block.
     * \param id ID of the building block.
     */
    TreeNode(shp<BuildingBlock> block,
             unsigned int id);

    /*! \brief Returns true if it is an outlet node.
     *
     * \return True if it is an outlet node.
     */
    bool isOutletNode() const;

    /*! \brief Returns true if it is an inlet node.
     *
     * \return True if it is an inlet node.
     */
    bool isInletNode() const;

    /*! \brief Returns true if it is an inlet or an outlet node.
     *
     * \return True if it is an inlet or an outlet node.
     */
    bool isExtremalNode() const;

    std::vector<GeometricFace> getOutlets() const;

    std::vector<shp<TreeNode>>              M_children;
    shp<BuildingBlock>                      M_block;
    unsigned int                            M_ID;
    unsigned int                            M_depth;
    unsigned int                            M_nChildren;
};

/// A tree containing multiple TreeNodes.
class TreeStructure
{
private:
    typedef shp<TreeNode>                  TreeNodePtr;
    typedef shp<BuildingBlock>             BuildingBlockPtr;
    typedef LifeV::VectorSmall<3>          Vector3D;
public:

    /*! \brief Constructor.
     *
     * \param verbose If true, output is pushed to standard output.
     */
    TreeStructure(bool verbose = false);

    /*! \brief Add child to specific node.
     *
     * \param baseID ID of the base node.
     * \param blockToAdd The block to add.
     * \param outletIndex The outlet to which the block must be added.
     * \return ID of the new node.
     */
    unsigned int addChild(unsigned int baseID,
                          BuildingBlockPtr blockToAdd,
                          int outletIndex = -1);

    /*! \brief Add child to specific node.
     *
     * \param baseID ID of the base node.
     * \param blockToAdd The block to add.
     * \param outletIndex The outlet to which the block must be added.
     */
    void addChild(unsigned int baseID,
                  TreeNodePtr nodeToAdd,
                  int outletIndex = -1);

    /*! \brief Set the root of the tree.
     *
     * \param Shared pointer to the building block of the head.
     */
    void setRoot(BuildingBlockPtr blockHead);

    /*! \brief Getter for the root node.
     *
     * \return Shared pointer to the root TreeNode.
     */
    TreeNodePtr getRoot();

    /*! \brief Get max ID of the nodes.
     *
     * \return The max ID.
     */
    unsigned int getMaxID();

    /*! \brief Returns true if the tree is empty.
     *
     * \return True if the tree is empty.
     */
    bool isEmpty();

    /*! \brief Returns maximum depth from root to any leaf.
     *
     * \return The desired length.
     */
    unsigned int depth();

    /*! \brief Apply geometric transformations to every node.
     *
     * \param deformMesh If true, the meshes in the nodes are deformed.
     */
    void traverseAndDeformGeometries(bool deformMesh = true);

    /// Call the same function on every inner BuildingBlock.
    void resetInletOutlets();

    /*! \brief Get set with all the mesh names.
     *
     * \return Standard set with all the mesh names.
     */
    std::set<std::string> getMeshListNames();

    /*! \brief Dump the tree to file.
     *
     * \param outdir The output directory.
     * \param meshdir The directory where meshes are stored.
     */
    void dump(std::string outdir,
              std::string meshdir);

    /*! \brief Read all the meshes.
     *
     * \param meshdir The directory where meshes are stored.
     */
    void readMeshes(std::string meshdir = "../geometries/");

    /*! \brief Create a random tree.
     *
     * \param blocksNumber The number of blocks.
     * \param comm The MPI Communicator.
     */
    void createRandom(unsigned int blocksNumber,
                      shp<Epetra_Comm> comm);

    /*! \brief Get a map with all the IDs and TreeNodes.
     *
     * \return Returns a map with key = ID and value = shared pointer to TreeNode.
     */
    std::map<unsigned int, TreeNodePtr> getNodesMap();

    /*! \brief Randomize geometric parameters of each building block.
     *
     * The parameters are sampled in an interval centered on the original value.
     *
     * \param bound The bounds around the original value.
     */
    void randomSampleAroundOriginalValue(const double& bound);

    /*! \brief Operator plus.
     *
     * Returns a TreeNode obtained as the union of the two operands. The trees
     * are joined if a suitable face is found.
     *
     * \param other The other TreeStructure.
     */
    TreeStructure& operator+(TreeStructure& other);

    void setGeometricParametersFromSample(const std::map<std::string, double>& sample);

    /*! \brief Search a block with a specific face.
     *
     * \param centerOfTheFace The center of the face to be looked for.
     * \param tol Tolerance to determine if two faces are the same.
     * \param outletIdx The index of the outlet.
     * \return The ID of the block; -1 if the block does not exist.
     */
    int findBlockWithFace(const Vector3D& centerOfTheFace,
                          const double& tol,
                          int& outletIdx);

    /// Recompute node IDs and depth of the tree.
    void resetNodesIDs();

private:

    TreeNodePtr                             M_root;
    unsigned int                            M_maxID;
    std::map<unsigned int, TreeNodePtr>     M_nodesMap;
    unsigned int                            M_depth;
    bool                                    M_verbose;
};

}  // namespace RedMA

#endif  // TREESTRUCTURE_HPP
