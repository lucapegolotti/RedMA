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

#include <BuildingBlock.hpp>
#include <memory>
#include <map>
#include <queue>

namespace RedMA
{

class TreeNode
{
public:
    TreeNode(std::shared_ptr<BuildingBlock> block, unsigned int id);

    std::vector<std::shared_ptr<TreeNode> > M_children;
    std::shared_ptr<BuildingBlock> M_block;
    unsigned int M_ID;
    unsigned int M_depth;
    unsigned int M_nChildren;
};

class TreeStructure
{
private:
    typedef std::shared_ptr<TreeNode> TreeNodePtr;
    typedef std::shared_ptr<BuildingBlock> BuildingBlockPtr;
public:
    TreeStructure(bool verbose = false);

    unsigned int addChild(unsigned int baseID,
                          BuildingBlockPtr blockToAdd, int outletIndex = -1);

    void setRoot(BuildingBlockPtr blockHead);

    unsigned int getMaxID();

    // void print();

    bool isEmpty();

    unsigned int depth();

    void traverseAndDeformGeometries();

    void dump(std::string outdir, std::string meshdir);

    void readMeshes(std::string meshdir = "../geometries/");

    void createRandom(unsigned int blocksNumber);

private:
    // std::vector<std::vector<std::string> > fillDepthVectors();

    TreeNodePtr M_root;
    unsigned int M_maxID;
    std::map<unsigned int, TreeNodePtr> M_nodesMap;
    unsigned int M_depth;
    bool M_verbose;
};

}  // namespace RedMA

#endif  // BUILDINGBLOCK_HPP
