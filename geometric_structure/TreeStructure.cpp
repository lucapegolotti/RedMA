#include <TreeStructure.hpp>

namespace ReMA
{

TreeNode::TreeNode(std::shared_ptr<BuildingBlock> block, unsigned int id) :
  M_block(block),
  M_ID(id)
{
}

TreeStructure::TreeStructure() :
  M_maxID(0)
{
}

unsigned int TreeStructure::getMaxID()
{
    return M_maxID;
}

unsigned int TreeStructure::addChild(unsigned int baseID, std::shared_ptr<BuildingBlock> blockToAdd)
{
    if (M_nodesMap.find(baseID) == M_nodesMap.end())
    {
        std::string errorMsg = "Node tree with ID = " + std::to_string(baseID) +
                               " does not exist!";
        throw Exception(errorMsg);
    }

    TreeNodePtr baseNode = M_nodesMap[baseID];
    if (baseNode->M_block->expectedNumberOfChildren() == baseNode->M_children.size())
    {
        std::string errorMsg = "Node tree with ID = " + std::to_string(baseNode->M_ID) +
                               " can not have other children!";
        throw Exception(errorMsg);
    }
    std::shared_ptr<TreeNode> newNode(new TreeNode(blockToAdd,M_maxID));
    baseNode->M_children.push_back(newNode);
    M_nodesMap[newNode->M_ID] = newNode;
    M_maxID++;
    return newNode->M_ID;
}

void TreeStructure::setRoot(std::shared_ptr<BuildingBlock> block)
{
    if (M_maxID != 0)
    {
        std::string errorMsg = "Trees can only contain one root!";
        throw Exception(errorMsg);
    }
    std::shared_ptr<TreeNode> newNode(new TreeNode(block,0));
    M_root = newNode;
    M_nodesMap[0] = M_root;
    M_maxID++;
}

bool TreeStructure::isEmpty()
{
    return (M_nodesMap.size() == 0);
}

}  // namespace ReMA
