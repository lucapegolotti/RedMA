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

// void
// TreeStructure::
// print()
// {
//     typedef std::vector<std::string> StringVector;
//     typedef std::vector<StringVector> StringVectorVector;
//     std::vector<StringVector> labels = fillDepthVectors();
//
//     for (StringVectorVector::iterator it = labels.begin();
//          it != labels.end(); it++)
//     {
//         for (StringVector::iterator jt = it->begin();
//              jt != it->end(); jt++)
//         {
//             printlog(MAGENTA, *jt + "\t");
//         }
//         printlog(WHITE, "\n");
//     }
// }

// std::vector<std::vector<std::string> >
// TreeStructure::
// fillDepthVectors()
// {
//     std::vector<std::vector<std::string> > returnVec(M_depth+1);
//     std::cout << "DEPTH = " << M_depth << std::endl;
//     std::queue<TreeNodePtr> nodesQueue;
//     nodesQueue.push(M_root);
//
//     returnVec[0].push_back(M_root->M_block->name());
//     while (nodesQueue.size() != 0)
//     {
//         TreeNodePtr curNode = nodesQueue.front();
//         nodesQueue.pop();
//         typedef std::vector<TreeNodePtr> TreeNodesVector;
//         TreeNodesVector& children = curNode->M_children;
//         unsigned int expectedChildren =
//                      curNode->M_block->expectedNumberOfChildren();
//         for (int i = 0; i < expectedChildren; i++)
//         {
//             TreeNodePtr curChild = children[i];
//             if (curChild != nullptr)
//             {
//                 std::cout << "curChild->M_depth = " << curChild->M_depth << std::endl;
//                 nodesQueue.push(curChild);
//                 returnVec[curChild->M_depth].push_back(curChild->M_block->name());
//             }
//             else
//             {
//                 std::cout << "curChild->M_depth+1 = " << curChild->M_depth+1 << std::endl;
//                 returnVec[curNode->M_depth].push_back("NULL");
//             }
//         }
//     };
//
//     return returnVec;
// }

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
        if (chosenClass == 0)
        {
            printlog(MAGENTA, "[TreeStructure] Initializing tube\n", M_verbose);
            newBlock.reset(new Tube(comm, M_verbose));
        }
        else
        {
            printlog(MAGENTA, std::string("[TreeStructure] Initializing ") +
                              " bifurcation symmetric\n", M_verbose);
            newBlock.reset(new BifurcationSymmetric(comm, M_verbose));
        }
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

}  // namespace RedMA
