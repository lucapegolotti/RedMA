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

#include <redma/geometry/BuildingBlock.hpp>
#include <redma/geometry/Tube.hpp>
#include <redma/geometry/BifurcationSymmetric.hpp>

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

class TreeNode
{
public:
    TreeNode(shp<BuildingBlock> block, unsigned int id);

    bool isOutletNode();

    bool isInletNode();

    std::vector<shp<TreeNode> > M_children;
    shp<BuildingBlock>          M_block;
    unsigned int                            M_ID;
    unsigned int                            M_depth;
    unsigned int                            M_nChildren;
};

class TreeStructure
{
private:
    typedef shp<TreeNode>      TreeNodePtr;
    typedef shp<BuildingBlock> BuildingBlockPtr;
    typedef LifeV::VectorSmall<3>          Vector3D;
public:
    TreeStructure(bool verbose = false);

    unsigned int addChild(unsigned int baseID,
                          BuildingBlockPtr blockToAdd, int outletIndex = -1);

    void addChild(unsigned int baseID, TreeNodePtr nodeToAdd, int outletIndex = -1);

    void setRoot(BuildingBlockPtr blockHead);

    TreeNodePtr getRoot();

    unsigned int getMaxID();

    bool isEmpty();

    unsigned int depth();

    void traverseAndDeformGeometries(bool deformMesh = true);

    void resetInletOutlets();

    std::set<std::string> getMeshListNames();

    void dump(std::string outdir, std::string meshdir);

    void readMeshes(std::string meshdir = "../geometries/");

    void createRandom(unsigned int blocksNumber, shp<Epetra_Comm> comm);

    std::map<unsigned int, TreeNodePtr> getNodesMap();

    void randomSampleAroundOriginalValue(const double& bound);

    TreeStructure& operator+(TreeStructure& other);

    int findBlockWithFace(const Vector3D& centerOfTheFace, const double& tol,
                          int& outletIdx);

    void resetNodesIDs();

private:

    TreeNodePtr                         M_root;
    unsigned int                        M_maxID;
    std::map<unsigned int, TreeNodePtr> M_nodesMap;
    unsigned int                        M_depth;
    bool                                M_verbose;
};

}  // namespace RedMA

#endif  // BUILDINGBLOCK_HPP
