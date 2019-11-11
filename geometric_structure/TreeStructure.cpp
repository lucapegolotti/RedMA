#include <TreeStructure.hpp>

namespace RedMA
{

TreeNode::
TreeNode(std::shared_ptr<BuildingBlock> block, unsigned int id) :
  M_block(block),
  M_ID(id),
  M_nChildren(0)
{
    M_children.resize(M_block->expectedNumberOfChildren());
    for (int i = 0; i < M_block->expectedNumberOfChildren(); i++)
        M_children[i] = nullptr;
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
addChild(unsigned int baseID, BuildingBlockPtr blockToAdd, int outletIndex)
{
    if (M_nodesMap.find(baseID) == M_nodesMap.end())
    {
        std::string errorMsg = "Node tree with ID = " + std::to_string(baseID) +
                               " does not exist!";
        throw Exception(errorMsg);
    }

    TreeNodePtr baseNode = M_nodesMap[baseID];
    // note: this ensures that 1) we don't add to many children when outlets
    // are not specifies in the xml file, 2) we don't overwrite a already added
    // child
    if (baseNode->M_block->expectedNumberOfChildren() == baseNode->M_nChildren)
    {
        std::string errorMsg = "Node tree with ID = " +
                               std::to_string(baseNode->M_ID) +
                               " can not have other children!";
        throw Exception(errorMsg);
    }
    if (outletIndex != -1 &&
        baseNode->M_block->expectedNumberOfChildren() <= outletIndex)
    {
        std::string maxChildrenStr =
                  std::to_string(baseNode->M_block->expectedNumberOfChildren());
        std::string errorMsg = "Node tree with ID = " +
                               std::to_string(baseNode->M_ID) +
                               " is not compatible with outletIndex " +
                               std::to_string(outletIndex) + " (max index = " +
                               maxChildrenStr + ")";
        throw Exception(errorMsg);
    }
    std::shared_ptr<TreeNode> newNode(new TreeNode(blockToAdd,M_maxID));
    if (outletIndex == -1)
        baseNode->M_children[baseNode->M_nChildren] = newNode;
    else
        baseNode->M_children[outletIndex] = newNode;

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
    baseNode->M_nChildren++;
    return newNode->M_ID;
}

void
TreeStructure::
addChild(unsigned int baseID, TreeNodePtr nodeToAdd, int outletIndex)
{
    if (M_nodesMap.find(baseID) == M_nodesMap.end())
    {
        std::string errorMsg = "Node tree with ID = " + std::to_string(baseID) +
                               " does not exist!";
        throw Exception(errorMsg);
    }

    TreeNodePtr baseNode = M_nodesMap[baseID];
    // note: this ensures that 1) we don't add too many children when outlets
    // are not specified in the xml file, 2) we don't overwrite a already added
    // child
    if (baseNode->M_block->expectedNumberOfChildren() == baseNode->M_nChildren)
    {
        std::string errorMsg = "Node tree with ID = " +
                               std::to_string(baseNode->M_ID) +
                               " can not have other children!";
        throw Exception(errorMsg);
    }
    if (outletIndex != -1 &&
        baseNode->M_block->expectedNumberOfChildren() <= outletIndex)
    {
        std::string maxChildrenStr =
                  std::to_string(baseNode->M_block->expectedNumberOfChildren());
        std::string errorMsg = "Node tree with ID = " +
                               std::to_string(baseNode->M_ID) +
                               " is not compatible with outletIndex " +
                               std::to_string(outletIndex) + " (max index = " +
                               maxChildrenStr + ")";
        throw Exception(errorMsg);
    }
    if (outletIndex == -1)
        baseNode->M_children[baseNode->M_nChildren] = nodeToAdd;
    else
        baseNode->M_children[outletIndex] = nodeToAdd;

    baseNode->M_nChildren++;
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
            TreeNodePtr curChild = children[i];

            if (curChild)
            {
                GeometricFace curFace = curNode->M_block->getOutlet(i);

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

void
TreeStructure::
createRandom(unsigned int blocksNumber, std::shared_ptr<Epetra_Comm> comm)
{
    srand(time(NULL));
    for (int i = 0; i < blocksNumber; i++)
    {
        const unsigned int numberClasses = 2;
        std::shared_ptr<BuildingBlock> newBlock;

        int chosenClass = rand() % numberClasses;
        std::string ref = "coarse";
        if (chosenClass == 0)
        {
            printlog(MAGENTA, "[TreeStructure] Initializing tube\n", M_verbose);
            newBlock.reset(new Tube(comm, ref, M_verbose));
        }
        else
        {
            printlog(MAGENTA, std::string("[TreeStructure] Initializing ") +
                              " bifurcation symmetric\n", M_verbose);
            newBlock.reset(new BifurcationSymmetric(comm, ref, M_verbose));
        }
        // else
        // {
        //     printlog(MAGENTA, std::string("[TreeStructure] Initializing ") +
        //                       " bifurcation asymmetric\n", M_verbose);
        //     newBlock.reset(new BifurcationAsymmetric(comm, ref, M_verbose));
        // }
        newBlock->setRandom();

        if (i == 0)
            setRoot(newBlock);
        else
        {
            int a = 0;
            for (auto it = M_nodesMap.begin(); it != M_nodesMap.end(); it++)
            {
                unsigned int id = it->first;
                if (it->second->M_nChildren <
                    it->second->M_block->expectedNumberOfChildren())
                {
                    addChild(id, newBlock);
                    break;
                }
            }
        }
    }
}

TreeStructure::TreeNodePtr
TreeStructure::
getRoot()
{
    return M_root;
}

std::map<unsigned int, TreeStructure::TreeNodePtr>
TreeStructure::
getNodesMap()
{
    return M_nodesMap;
}

int
TreeStructure::
findBlockWithFace(const Vector3D& centerOfTheFace, const double& tol, int& outletIdx)
{
    std::queue<TreeNodePtr> nodesQueue;
    nodesQueue.push(M_root);
    while (nodesQueue.size() != 0)
    {
        TreeNodePtr curNode = nodesQueue.front();
        nodesQueue.pop();

        if (curNode->M_depth == 0)
        {
            if ((curNode->M_block->getInlet().M_center-centerOfTheFace).norm() < tol)
                return curNode->M_ID;
        }

        if (curNode->M_nChildren == 0)
        {
            std::vector<GeometricFace> outlets = curNode->M_block->getOutlets();

            unsigned int count = 0;
            for (auto it = outlets.begin(); it != outlets.end(); it++)
            {
                if ((it->M_center - centerOfTheFace).norm() < tol)
                {
                    outletIdx = count;
                    return curNode->M_ID;
                }
                count++;
            }
        }


        typedef std::vector<TreeNodePtr> TreeNodesVector;
        TreeNodesVector& children = curNode->M_children;
        unsigned int expectedChildren =
                     curNode->M_block->expectedNumberOfChildren();
        for (int i = 0; i < expectedChildren; i++)
        {
            TreeNodePtr curChild = children[i];

            if (curChild)
                nodesQueue.push(curChild);
        }
    }

    // return -1 if face does not belong
    return -1;
}

TreeStructure&
TreeStructure::
operator+=(TreeStructure& other)
{
    if (M_root == nullptr)
    {
        M_root = other.M_root;
        return *this;
    }

    const double tol = 1e-12;
    int status;
    int outletIndex = -1;
    // first we see if other must be attached to this
    status = findBlockWithFace(other.M_root->M_block->getInlet().M_center,
                               tol, outletIndex);

    if (status > 0)
    {
        addChild(status, other.M_root, outletIndex);
        return *this;
    }
    else
    {
        status = other.findBlockWithFace(M_root->M_block->getInlet().M_center,
                                         tol, outletIndex);
        if (status <= 0)
            throw new Exception("Operator += of treeStructure: no suitable face!");

        other.addChild(status, M_root, outletIndex);

        other.resetNodesIDs();
        return other;
    }
}

void
TreeStructure::
resetNodesIDs()
{
    std::queue<TreeNodePtr> nodesQueue;
    nodesQueue.push(M_root);
    unsigned int count = 0;
    M_root->M_depth = 0;
    M_depth = 0;
    while (nodesQueue.size() != 0)
    {
        TreeNodePtr curNode = nodesQueue.front();
        nodesQueue.pop();

        curNode->M_ID = count;

        typedef std::vector<TreeNodePtr> TreeNodesVector;
        TreeNodesVector& children = curNode->M_children;
        unsigned int expectedChildren =
                     curNode->M_block->expectedNumberOfChildren();
        for (int i = 0; i < expectedChildren; i++)
        {
            TreeNodePtr curChild = children[i];
            if (curChild)
            {
                curChild->M_depth = curNode->M_depth+1;
                M_depth = M_depth > curChild->M_depth ? M_depth :
                                                        curChild->M_depth;
                nodesQueue.push(curChild);
            }
        }
        count++;
    }
}


}  // namespace RedMA
