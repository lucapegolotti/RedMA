#include <TreeStructure.hpp>

namespace ReMA
{

TreeNode::
TreeNode(std::shared_ptr<BuildingBlock> block, unsigned int id) :
  M_block(block),
  M_ID(id)
{
}

TreeStructure::
TreeStructure() :
  M_maxID(0)
{
}

unsigned int
TreeStructure::
getMaxID()
{
    return M_maxID;
}

unsigned int
TreeStructure::
addChild(unsigned int baseID, BuildingBlockPtr blockToAdd)
{
    if (M_nodesMap.find(baseID) == M_nodesMap.end())
    {
        std::string errorMsg = "Node tree with ID = " + std::to_string(baseID) +
                               " does not exist!";
        throw Exception(errorMsg);
    }

    TreeNodePtr baseNode = M_nodesMap[baseID];
    if (baseNode->M_block->expectedNumberOfChildren() ==
        baseNode->M_children.size())
    {
        std::string errorMsg = "Node tree with ID = " +
                               std::to_string(baseNode->M_ID) +
                               " can not have other children!";
        throw Exception(errorMsg);
    }
    std::shared_ptr<TreeNode> newNode(new TreeNode(blockToAdd,M_maxID));
    baseNode->M_children.push_back(newNode);
    newNode->M_depth = baseNode->M_depth + 1;
    M_nodesMap[newNode->M_ID] = newNode;
    M_maxID++;
    M_depth = newNode->M_depth > M_depth ? newNode->M_depth : M_depth;
    return newNode->M_ID;
}

void
TreeStructure::
setRoot(BuildingBlockPtr block)
{
    if (M_maxID != 0)
    {
        std::string errorMsg = "Trees can only contain one root!";
        throw Exception(errorMsg);
    }
    std::shared_ptr<TreeNode> newNode(new TreeNode(block,0));
    M_root = newNode;
    M_root->M_depth = 0;
    M_nodesMap[0] = M_root;
    M_maxID++;
    M_depth = 0;
}

bool
TreeStructure::
isEmpty()
{
    return (M_nodesMap.size() == 0);
}

void
TreeStructure::
print()
{
    typedef std::vector<std::string> StringVector;
    typedef std::vector<StringVector> StringVectorVector;
    std::vector<StringVector> labels = fillDepthVectors();

    for (StringVectorVector::iterator it = labels.begin();
         it != labels.end(); it++)
    {
        for (StringVector::iterator jt = it->begin();
             jt != it->end(); jt++)
        {
            printlog(MAGENTA, *jt + "\t");
        }
        printlog(WHITE, "\n");
    }
}

std::vector<std::vector<std::string> >
TreeStructure::
fillDepthVectors()
{
    std::vector<std::vector<std::string> > returnVec(M_depth+1);

    TreeNodePtr curNode = M_root;
    std::queue<TreeNodePtr > nodesQueue;

    returnVec[0].push_back(M_root->M_block->name());
    do
    {
        typedef std::vector<TreeNodePtr> TreeNodesVector;
        TreeNodesVector& children = curNode->M_children;
        unsigned int expectedChildren =
                     curNode->M_block->expectedNumberOfChildren();
        for (int i = 0; i < expectedChildren; i++)
        {
            if (i < children.size())
            {
                TreeNodePtr curChild = children[i];
                nodesQueue.push(curChild);
                returnVec[curChild->M_depth].push_back(curChild->M_block->name());
            }
            else
            {
                returnVec[curNode->M_depth+1].push_back("NULL");
            }
        }
        curNode = nodesQueue.front();
        nodesQueue.pop();
    } while (nodesQueue.size() != 0);

    return returnVec;
}

unsigned int
TreeStructure::
depth()
{
    return M_depth;
}

}  // namespace ReMA
