#include <TreeStructure.hpp>

namespace RedMA
{

TreeNode::
TreeNode(std::shared_ptr<BuildingBlock> block, unsigned int id) :
  M_block(block),
  M_ID(id)
{
}

TreeStructure::
TreeStructure(bool verbose) :
  M_maxID(0),
  M_verbose(verbose)
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
    if (M_nodesMap.find(newNode->M_ID) != M_nodesMap.end())
    {
        std::string errorMsg = "Node tree with ID = " +
                               std::to_string(newNode->M_ID) +
                               " already exists!";
        throw Exception(errorMsg);
    }
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
    std::queue<TreeNodePtr> nodesQueue;
    nodesQueue.push(M_root);

    returnVec[0].push_back(M_root->M_block->name());
    while (nodesQueue.size() != 0)
    {
        TreeNodePtr curNode = nodesQueue.front();
        nodesQueue.pop();
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
    };

    return returnVec;
}

void
TreeStructure::
traverseAndDeformGeometries()
{
    std::queue<TreeNodePtr> nodesQueue;
    nodesQueue.push(M_root);
    while(nodesQueue.size() != 0)
    {
        TreeNodePtr curNode = nodesQueue.front();
        nodesQueue.pop();
        curNode->M_block->applyGlobalTransformation();
        typedef std::vector<TreeNodePtr> TreeNodesVector;
        TreeNodesVector& children = curNode->M_children;
        unsigned int expectedChildren =
                     curNode->M_block->expectedNumberOfChildren();
        for (int i = 0; i < expectedChildren; i++)
        {
            if (i < children.size())
            {
                GeometricFace curFace = curNode->M_block->getOutlet(i);

                TreeNodePtr curChild = children[i];
                curChild->M_block->mapChildInletToParentOutlet(curFace);
                nodesQueue.push(curChild);
            }
        }
    }
}

unsigned int
TreeStructure::
depth()
{
    return M_depth;
}

void
TreeStructure::
dump(std::string outdir, std::string meshdir)
{
    typedef std::map<unsigned int, TreeNodePtr> NodesMap;

    for (NodesMap::iterator it = M_nodesMap.begin();
         it != M_nodesMap.end(); it++)
    {
        std::string nodeName = (*it).second->M_block->name() +
                               std::to_string((*it).second->M_ID);
        (*it).second->M_block->dumpMesh(outdir, meshdir, nodeName);
    }
}

void
TreeStructure::
readMeshes(std::string meshdir)
{
    typedef std::map<unsigned int, TreeNodePtr> NodesMap;

    for (NodesMap::iterator it = M_nodesMap.begin();
         it != M_nodesMap.end(); it++)
    {
        (*it).second->M_block->readMesh(meshdir);
    }
}

}  // namespace RedMA
