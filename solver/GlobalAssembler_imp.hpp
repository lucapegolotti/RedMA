// implementation of template class

namespace RedMA
{

template <class AssemblerType>
GlobalAssembler<AssemblerType>::
GlobalAssembler(const GetPot& datafile, commPtr_Type comm, bool verbose) :
  M_datafile(datafile),
  M_comm(comm),
  M_verbose(verbose)
{
    M_globalMap.reset(new MapEpetra());
}

template <class AssemblerType>
void
GlobalAssembler<AssemblerType>::
buildPrimalStructures(TreeStructure& tree)
{
    typedef std::map<unsigned int, TreeNodePtr>     NodesMap;

    NodesMap nodesMap = tree.getNodesMap();

    for (NodesMap::iterator it = nodesMap.begin(); it != nodesMap.end(); it++)
    {
        AssemblerTypePtr newAssembler(new AssemblerType(M_datafile, M_comm,
                                                        it->second, M_verbose));
        newAssembler->setup();

        newAssembler->addPrimalMaps(M_globalMap, M_dimensionsVector);
        M_assemblersVector.push_back(std::make_pair(it->first, newAssembler));
        M_assemblersMap[it->first] = newAssembler;
    }
}

template <class AssemblerType>
void
GlobalAssembler<AssemblerType>::
setup(TreeStructure& tree)
{
    buildPrimalStructures(tree);
    buildDualStructures(tree);

    M_offsets.push_back(0);

    unsigned int offset = 0;
    for (std::vector<unsigned int>::iterator it = M_dimensionsVector.begin();
         it != M_dimensionsVector.end(); it++)
    {
        unsigned int val = *it;
        M_offsets.push_back(offset + val);
        offset += val;
    }
}

template <class AssemblerType>
void
GlobalAssembler<AssemblerType>::
buildDualStructures(TreeStructure& tree)
{
    typedef std::map<unsigned int, TreeNodePtr>                     NodesMap;
    typedef std::vector<TreeNodePtr>                                NodesVector;
    typedef std::vector<std::pair<unsigned int, AssemblerTypePtr> > AssemblersVector;

    NodesMap nodesMap = tree.getNodesMap();

    // construct a vector of pairs for the indices of the domains sharing
    // interfaces.
    // we identify the interface with the index in the vector
    unsigned int interfaceCount = 0;
    for (NodesMap::iterator it = nodesMap.begin(); it != nodesMap.end(); it++)
    {
        NodesVector children = it->second->M_children;

        unsigned int countOutlet = 0;
        unsigned int myID = it->second->M_ID;
        AssemblerTypePtr fatherAssembler;

        fatherAssembler = M_assemblersMap[myID];
        for (NodesVector::iterator itVector = children.begin();
             itVector != children.end(); itVector++)
        {
            if (*itVector)
            {
                unsigned int otherID = (*itVector)->M_ID;
                M_interfaces.push_back(std::make_pair(myID, otherID));

                AssemblerTypePtr childAssembler;
                childAssembler = M_assemblersMap[otherID];
                // within this function, we also add to the global maps
                // the newly created map for the lagrange multiplier
                fatherAssembler->assembleCouplingMatrices(*childAssembler,
                                                        countOutlet,
                                                        interfaceCount,
                                                        M_globalMap,
                                                        M_dimensionsVector);
            }
            countOutlet++;
            interfaceCount++;
        }
    }
}

template <class AssemblerType>
typename GlobalAssembler<AssemblerType>::MapEpetraPtr
GlobalAssembler<AssemblerType>::
getGlobalMap() const
{
    return M_globalMap;
}

template <class AssemblerType>
typename GlobalAssembler<AssemblerType>::MatrixPtr
GlobalAssembler<AssemblerType>::
getGlobalMass() const
{
    return M_massMatrix;
}

template <class AssemblerType>
typename GlobalAssembler<AssemblerType>::MatrixPtr
GlobalAssembler<AssemblerType>::
getJacobianF(bool addCoupling, double* diagonalCoefficient)
{
    MatrixPtr jacobian;
    fillGlobalMatrix(jacobian, addCoupling, &AssemblerType::getJacobian,
                     diagonalCoefficient);
    return jacobian;
}

template <class AssemblerType>
typename GlobalAssembler<AssemblerType>::VectorPtr
GlobalAssembler<AssemblerType>::
computeF()
{
    VectorPtr f(new Vector(*M_globalMap));
    fillGlobalVector(f, &AssemblerType::computeF);
    return f;
}

template <class AssemblerType>
typename GlobalAssembler<AssemblerType>::VectorPtr
GlobalAssembler<AssemblerType>::
computeFder()
{
    VectorPtr fder(new Vector(*M_globalMap));
    fillGlobalVector(fder, &AssemblerType::computeFder);
    return fder;
}

template <class AssemblerType>
typename GlobalAssembler<AssemblerType>::MatrixPtr
GlobalAssembler<AssemblerType>::
assembleGlobalMass(bool addCoupling, double* diagonalCoefficient)
{
    fillGlobalMatrix(M_massMatrix, addCoupling, &AssemblerType::getMassMatrix,
                     diagonalCoefficient);
    return M_massMatrix;
}

template<class AssemblerType>
template<typename FunctionType>
void
GlobalAssembler<AssemblerType>::
fillGlobalMatrix(MatrixPtr& matrixToFill, bool addCoupling,
                 FunctionType getMatrixMethod, double* diagonalCoefficient)
{
    using namespace LifeV::MatrixEpetraStructuredUtility;

    typedef std::vector<std::pair<unsigned int, AssemblerTypePtr> >
                AssemblersVector;

    typedef LifeV::MatrixEpetraStructuredView<double> MatrixView;

    matrixToFill.reset(new Matrix(*M_globalMap));
    matrixToFill->zero();

    LifeV::MatrixBlockStructure structure;
    structure.setBlockStructure(M_dimensionsVector,
                                M_dimensionsVector);

    // we start with the primal blocks
    unsigned int countBlocks = 0;
    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        unsigned int blockIndex = it->first;
        unsigned int numberBlocks = it->second->numberOfBlocks();
        for (int i = 0; i < numberBlocks; i++)
        {
            for (int j = 0; j < numberBlocks; j++)
            {
                AssemblerType& curAssembler = *it->second;
                MatrixPtr localMatrix = (curAssembler.*getMatrixMethod)(i, j);
                if (diagonalCoefficient)
                    curAssembler.applyBCsMatrix(localMatrix, *diagonalCoefficient,
                                                i, j);
                if (localMatrix)
                {
                    std::shared_ptr<MatrixView> blockGlobalView;
                    blockGlobalView = createBlockView(matrixToFill, structure,
                                                      i + countBlocks,
                                                      j + countBlocks);

                    LifeV::MatrixBlockStructure blockStructure;
                    std::vector<unsigned int> rows(1), cols(1);
                    rows[0] = M_dimensionsVector[i + countBlocks];
                    cols[0] = M_dimensionsVector[j + countBlocks];
                    blockStructure.setBlockStructure(rows, cols);

                    std::shared_ptr<MatrixView> blockLocalView;
                    blockLocalView =
                        createBlockView(localMatrix, blockStructure, 0, 0);
                    copyBlock(blockLocalView, blockGlobalView);
                }
            }
        }
        countBlocks += numberBlocks;
    }

    // if required add the coupling blocks
    if (addCoupling)
    {
        typedef std::vector<std::pair<unsigned int, unsigned int> >
                InterfacesVector;

        unsigned int offset = 0;
        for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
             it != M_assemblersVector.end(); it++)
             offset += it->second->numberOfBlocks();

        for (InterfacesVector::iterator it = M_interfaces.begin();
             it != M_interfaces.end(); it++)
        {
            unsigned int indices[2] = {it->first, it->second};

            for (int i = 0; i < 2; i++)
            {
                AssemblerType& curAssembler = *M_assemblersMap[indices[i]];
                MatrixPtr Qt = curAssembler.getQT(indices[(i+1) % 2]);
                unsigned int blockCoupling = curAssembler.getIndexCoupling();

                curAssembler.applyBCsMatrix(Qt, 0.0,
                                            blockCoupling, blockCoupling);

                unsigned int blockIndex = blockCoupling +
                                indices[i] * curAssembler.numberOfBlocks();

                if (Qt)
                {
                    std::shared_ptr<MatrixView> blockGlobalView;
                    blockGlobalView = createBlockView(matrixToFill, structure,
                                                      blockIndex, offset);

                    LifeV::MatrixBlockStructure blockStructure;
                    std::vector<unsigned int> rows(1), cols(1);
                    rows[0] = M_dimensionsVector[blockIndex];
                    cols[0] = M_dimensionsVector[offset];
                    blockStructure.setBlockStructure(rows, cols);

                    std::shared_ptr<MatrixView> blockLocalView;
                    blockLocalView = createBlockView(Qt, blockStructure, 0, 0);
                    copyBlock(blockLocalView, blockGlobalView);
                }

                MatrixPtr Q = curAssembler.getQ(indices[(i+1) % 2]);

                if (Q)
                {
                    std::shared_ptr<MatrixView> blockGlobalView;
                    blockGlobalView = createBlockView(matrixToFill, structure,
                                                      offset, blockIndex);

                    LifeV::MatrixBlockStructure blockStructure;
                    std::vector<unsigned int> rows(1), cols(1);
                    cols[0] = M_dimensionsVector[blockIndex];
                    rows[0] = M_dimensionsVector[offset];
                    blockStructure.setBlockStructure(rows, cols);

                    std::shared_ptr<MatrixView> blockLocalView;
                    blockLocalView = createBlockView(Q, blockStructure, 0, 0);
                    copyBlock(blockLocalView, blockGlobalView);
                }
            }
            offset++;
        }
    }

    matrixToFill->globalAssemble();
}

template<class AssemblerType>
template<typename FunctionType>
void
GlobalAssembler<AssemblerType>::
fillGlobalVector(VectorPtr& vectorToFill, FunctionType getVectorMethod)
{
    using namespace LifeV::MatrixEpetraStructuredUtility;

    typedef std::vector<std::pair<unsigned int, AssemblerTypePtr> >
                AssemblersVector;

    typedef std::vector<MapEpetraPtr>                    MapVector;

    vectorToFill.reset(new Vector(*M_globalMap));
    vectorToFill->zero();

    unsigned int nAssemblers = M_assemblersVector.size();
    unsigned int offset = 0;
    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        std::vector<VectorPtr> localSolutions;
        MapVector maps = it->second->getPrimalMapVector();
        AssemblerType& curAssembler = *it->second;
        std::vector<VectorPtr> localVectors = (curAssembler.*getVectorMethod)();
        unsigned int index = 0;
        for (MapVector::iterator itmap = maps.begin();
             itmap != maps.end(); itmap++)
        {
            LifeV::MapEpetra& curLocalMap = **itmap;
            // we fill only if the vector corresponding to map i exists!
            if (localVectors[index])
                vectorToFill->subset(*localVectors[index], curLocalMap, 0,
                                     offset);
            offset += curLocalMap.mapSize();
            index++;
        }

        // deal with the dual part. This is trickier because more assemblers
        // contribute to the same subsets of vectofill
        index = 0;
        maps = it->second->getDualMapVector();
        std::vector<unsigned int> indices = it->second->getInterfacesIndices();
        for (MapVector::iterator itmap = maps.begin();
             itmap != maps.end(); itmap++)
        {
            LifeV::MapEpetra& curLocalMap = **itmap;
            VectorPtr aux(new Vector(curLocalMap));
            aux->subset(*vectorToFill, curLocalMap,
                        M_offsets[nAssemblers + indices[index]], 0);
            *aux += *localVectors[index];
            vectorToFill->subset(*aux, curLocalMap, 0,
                                 M_offsets[nAssemblers + indices[index]]);
            index++;
        }
    }
}

template<class AssemblerType>
void
GlobalAssembler<AssemblerType>::
setTimeAndPrevSolution(const double& time, VectorPtr solution)
{
    typedef std::pair<unsigned int, AssemblerTypePtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;
    typedef std::shared_ptr<LifeV::MapEpetra>            MapEpetraPtr;
    typedef std::vector<MapEpetraPtr>                    MapVector;

    unsigned int offset = 0;
    unsigned int nAssemblers = M_assemblersVector.size();
    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        std::vector<VectorPtr> localSolutions;
        // first handle the primal solutions
        MapVector maps = it->second->getPrimalMapVector();
        for (MapVector::iterator itmap = maps.begin();
             itmap != maps.end(); itmap++)
        {
            LifeV::MapEpetra& curLocalMap = **itmap;
            VectorPtr subSolution;
            subSolution.reset(new Vector(curLocalMap));
            subSolution->zero();
            subSolution->subset(*solution, curLocalMap, offset, 0);
            localSolutions.push_back(subSolution);
            offset += curLocalMap.mapSize();
        }

        maps = it->second->getDualMapVector();
        std::vector<unsigned int> indices = it->second->getInterfacesIndices();
        unsigned int in = 0;
        for (MapVector::iterator itmap = maps.begin();
             itmap != maps.end(); itmap++)
        {
            LifeV::MapEpetra& curLocalMap = **itmap;
            VectorPtr subSolution;
            subSolution.reset(new Vector(curLocalMap));
            subSolution->zero();
            subSolution->subset(*solution, curLocalMap,
                                M_offsets[nAssemblers + indices[in]], 0);
            localSolutions.push_back(subSolution);
            offset += curLocalMap.mapSize();
            in++;
        }
        it->second->setTimeAndPrevSolution(time, localSolutions);
    }
}

template<class AssemblerType>
void
GlobalAssembler<AssemblerType>::
applyBCsRhsRosenbrock(VectorPtr rhs, VectorPtr utilde,
                      const double& time, const double& dt,
                      const double& alphai, const double& gammai)
{
    typedef std::pair<unsigned int, AssemblerTypePtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;
    typedef std::shared_ptr<LifeV::MapEpetra>            MapEpetraPtr;
    typedef std::vector<MapEpetraPtr>                    MapVector;

    unsigned int offset = 0;
    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        std::vector<VectorPtr> rhss;
        std::vector<VectorPtr> utildes;
        MapVector maps = it->second->getPrimalMapVector();
        unsigned int suboffset = 0;
        for (MapVector::iterator itmap = maps.begin();
             itmap != maps.end(); itmap++)
        {
            LifeV::MapEpetra& curLocalMap = **itmap;
            VectorPtr subRhs;
            VectorPtr subUtilde;
            subRhs.reset(new Vector(curLocalMap));
            subUtilde.reset(new Vector(curLocalMap));
            subRhs->zero();
            subUtilde->zero();
            subRhs->subset(*rhs, curLocalMap, offset + suboffset, 0);
            subUtilde->subset(*utilde, curLocalMap, offset + suboffset, 0);
            rhss.push_back(subRhs);
            utildes.push_back(subUtilde);
            suboffset += curLocalMap.mapSize();
        }
        // apply bcs
        it->second->applyBCsRhsRosenbrock(rhss, utildes, time, dt,
                                          alphai, gammai);

        suboffset = 0;
        unsigned int count = 0;
        // copy back to global vectors
        for (MapVector::iterator itmap = maps.begin();
             itmap != maps.end(); itmap++)
        {
          LifeV::MapEpetra& curLocalMap = **itmap;
          rhs->subset(*rhss[count], curLocalMap, 0, offset + suboffset);
          suboffset += curLocalMap.mapSize();
          count++;
        }
        offset += suboffset;
    }
}

template<class AssemblerType>
void
GlobalAssembler<AssemblerType>::
setMaxVelocityLawInflow(std::function<double(double)> maxLaw)
{
    typedef std::pair<unsigned int, AssemblerTypePtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;

    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        it->second->setMaxVelocityLawInflow(maxLaw);
    }
}

template<class AssemblerType>
void
GlobalAssembler<AssemblerType>::
setMaxVelocityDtLawInflow(std::function<double(double)> maxLawDt)
{
    typedef std::pair<unsigned int, AssemblerTypePtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;

    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        it->second->setMaxVelocityDtLawInflow(maxLawDt);
    }
}

template<class AssemblerType>
void
GlobalAssembler<AssemblerType>::
exportSolutions(const double& time, VectorPtr solution)
{
    typedef std::pair<unsigned int, AssemblerTypePtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;
    typedef std::shared_ptr<LifeV::MapEpetra>            MapEpetraPtr;
    typedef std::vector<MapEpetraPtr>                    MapVector;

    std::string solutionsDir;
    solutionsDir = M_datafile("exporter/outdirectory", "solution");

    unsigned int offset = 0;
    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        std::vector<VectorPtr> localSolutions;
        MapVector maps = it->second->getPrimalMapVector();
        for (MapVector::iterator itmap = maps.begin();
             itmap != maps.end(); itmap++)
        {
            LifeV::MapEpetra& curLocalMap = **itmap;
            VectorPtr subSolution;
            subSolution.reset(new Vector(curLocalMap));
            subSolution->zero();
            subSolution->subset(*solution, curLocalMap, offset, 0);
            localSolutions.push_back(subSolution);
            offset += curLocalMap.mapSize();
        }
        it->second->exportSolutions(time, localSolutions);
    }
}


}  // namespace RedMA
