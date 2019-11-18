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

        newAssembler->addPrimalMaps(M_globalMap, M_maps, M_dimensionsVector);
        M_assemblersVector.push_back(std::make_pair(it->first, newAssembler));
        M_assemblersMap[it->first] = newAssembler;
    }
    unsigned int count = 0;
    for (auto it = M_dimensionsVector.begin(); it != M_dimensionsVector.end(); it++)
        count += *it;
    std::string msg = "[GlobalAssembler] number of primal dofs = ";
    msg += std::to_string(count);
    msg += "\n";
    printlog(MAGENTA, msg, M_verbose);
}

template <class AssemblerType>
void
GlobalAssembler<AssemblerType>::
setup(TreeStructure& tree)
{
    typedef std::vector<std::pair<unsigned int, AssemblerTypePtr> > AssemblersVector;

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

    M_nPrimalBlocks = 0;
    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        M_nPrimalBlocks += it->second->numberOfBlocks();
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
                                                          M_maps,
                                                          M_dimensionsVector);
                interfaceCount++;
            }
            countOutlet++;
        }
    }
    unsigned int count = 0;
    for (auto it = M_dimensionsVector.begin(); it != M_dimensionsVector.end(); it++)
        count += *it;
    std::string msg = "[GlobalAssembler] number of primal + dual dofs = ";
    msg += std::to_string(count);
    msg += "\n";
    printlog(MAGENTA, msg, M_verbose);
}

template <class AssemblerType>
typename GlobalAssembler<AssemblerType>::MapEpetraPtr
GlobalAssembler<AssemblerType>::
getGlobalMap() const
{
    return M_globalMap;
}

template <class AssemblerType>
GlobalBlockMatrix
GlobalAssembler<AssemblerType>::
getGlobalMass()
{
    GlobalBlockMatrix updatedMassMatrix(M_massMatrix.getNumberRows(),
                                        M_massMatrix.getNumberCols());
    double diagCoefficient = 0.0;
    fillGlobalMatrix(updatedMassMatrix, false, &AssemblerType::getUpdateMass,
                     &diagCoefficient);
    updatedMassMatrix.add(M_massMatrix);
    return updatedMassMatrix;
    // return M_massMatrix;
}
template <class AssemblerType>
typename GlobalAssembler<AssemblerType>::VectorPtr
GlobalAssembler<AssemblerType>::
getInitialCondition()
{
    VectorPtr f(new Vector(*M_globalMap));
    fillGlobalVector(f, &AssemblerType::initialCondition);
    return f;
}

template <class AssemblerType>
GlobalBlockMatrix
GlobalAssembler<AssemblerType>::
getGlobalMassJac()
{
    GlobalBlockMatrix updatedMassMatrix(M_massMatrix.getNumberRows(),
                                        M_massMatrix.getNumberCols());
    double diagCoefficient = 0.0;
    fillGlobalMatrix(updatedMassMatrix, false, &AssemblerType::getUpdateMassJac,
                     &diagCoefficient);
    updatedMassMatrix.add(M_massMatrix);
    // return updatedMassMatrix;
    return updatedMassMatrix;
}

template <class AssemblerType>
GlobalBlockMatrix
GlobalAssembler<AssemblerType>::
getGlobalMassJacVelocity()
{
    GlobalBlockMatrix updatedMassMatrix(M_massMatrix.getNumberRows(),
                                        M_massMatrix.getNumberCols());
    double diagCoefficient = 0.0;
    fillGlobalMatrix(updatedMassMatrix, false, &AssemblerType::getUpdateMassJacVelocity,
                     &diagCoefficient);
    return updatedMassMatrix;
}

template <class AssemblerType>
GlobalBlockMatrix
GlobalAssembler<AssemblerType>::
getJacobianF(bool addCoupling, double* diagonalCoefficient)
{
    GlobalBlockMatrix jacobian;
    fillGlobalMatrix(jacobian, addCoupling, &AssemblerType::getJacobian,
                     diagonalCoefficient);
    return jacobian;
}

template <class AssemblerType>
GlobalBlockMatrix
GlobalAssembler<AssemblerType>::
getJacobianFprec(bool addCoupling, double* diagonalCoefficient)
{
    GlobalBlockMatrix jacobian;
    fillGlobalMatrix(jacobian, addCoupling, &AssemblerType::getJacobianPrec,
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
GlobalBlockMatrix
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
fillGlobalMatrix(GlobalBlockMatrix& matrixToFill, bool addCoupling,
                 FunctionType getMatrixMethod, double* diagonalCoefficient)
{
    typedef std::vector<std::pair<unsigned int, AssemblerTypePtr> >
                AssemblersVector;

    unsigned int totalNumberBlocks = M_nPrimalBlocks + M_interfaces.size();
    matrixToFill.resize(totalNumberBlocks, totalNumberBlocks);
    matrixToFill.setMaps(M_maps, M_maps);
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
                // Attention: this does not work if number of blocks is not constant
                // over all the domains
                matrixToFill.copyBlock(countBlocks + i,
                                       countBlocks + j,
                                       localMatrix);
            }
        }
        countBlocks += numberBlocks;
    }

    // if required add the coupling blocks
    if (addCoupling)
    {
        typedef std::vector<std::pair<unsigned int, unsigned int> >
                InterfacesVector;

        unsigned int offset = M_nPrimalBlocks;

        for (InterfacesVector::iterator it = M_interfaces.begin();
             it != M_interfaces.end(); it++)
        {
            unsigned int indices[2] = {it->first, it->second};

            for (int i = 0; i < 2; i++)
            {
                AssemblerType& curAssembler = *M_assemblersMap[indices[i]];
                MatrixPtr Qt = curAssembler.getQT(indices[(i+1) % 2]);
                unsigned int blockCoupling = curAssembler.getIndexCoupling();

                unsigned int blockIndex = blockCoupling +
                                indices[i] * curAssembler.numberOfBlocks();

                curAssembler.applyBCsMatrix(Qt, 0.0,
                                            blockCoupling, blockCoupling);
                matrixToFill.copyBlock(blockIndex, offset, Qt);

                MatrixPtr Q = curAssembler.getQ(indices[(i+1) % 2]);
                matrixToFill.copyBlock(offset, blockIndex, Q);
            }
            offset++;
        }
    }
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

    unsigned int offset = 0;
    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        std::vector<VectorPtr> localSolutions;
        MapVector maps = it->second->getPrimalMapVector();

        AssemblerType& curAssembler = *it->second;
        std::vector<VectorPtr> localVectors;
        localVectors = (curAssembler.*getVectorMethod)();

        unsigned int index = 0;

        for (MapVector::iterator itmap = maps.begin();
             itmap != maps.end(); itmap++)
        {
            LifeV::MapEpetra& curLocalMap = **itmap;
            // we fill only if the vector corresponding to map i exists!
            if (localVectors[index])
            {
                vectorToFill->subset(*localVectors[index], curLocalMap, 0,
                                     offset);
            }
            offset += curLocalMap.mapSize();
            index++;
        }

        // deal with the dual part. This is trickier because more assemblers
        // contribute to the same subsets of vectorfill
        index = 0;
        maps = it->second->getDualMapVector();
        std::vector<unsigned int> indices = it->second->getInterfacesIndices();
        for (MapVector::iterator itmap = maps.begin();
             itmap != maps.end(); itmap++)
        {
            LifeV::MapEpetra& curLocalMap = **itmap;
            VectorPtr aux(new Vector(curLocalMap));
            aux->subset(*vectorToFill, curLocalMap,
                        M_offsets[M_nPrimalBlocks + indices[index]], 0);
            *aux += 0;
            localVectors[index + it->second->numberOfBlocks()]->zero();
            vectorToFill->subset(*aux, curLocalMap, 0,
                                  M_offsets[M_nPrimalBlocks + indices[index]]);
            index++;
        }
    }
}

template<class AssemblerType>
void
GlobalAssembler<AssemblerType>::
setTimeAndPrevSolution(const double& time, VectorPtr solution, bool doAssembly)
{
    typedef std::pair<unsigned int, AssemblerTypePtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;
    typedef std::shared_ptr<LifeV::MapEpetra>            MapEpetraPtr;
    typedef std::vector<MapEpetraPtr>                    MapVector;

    unsigned int offsetPrimal = 0;
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
            subSolution->subset(*solution, curLocalMap, offsetPrimal, 0);
            localSolutions.push_back(subSolution);
            offsetPrimal += curLocalMap.mapSize();
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
                                 M_offsets[M_nPrimalBlocks + indices[in]], 0);
            localSolutions.push_back(subSolution);
            in++;
        }
        it->second->setTimeAndPrevSolution(time, localSolutions, doAssembly);
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
template<typename FunctionType>
void
GlobalAssembler<AssemblerType>::
applyBCsVector(VectorPtr rhs, const double& coeff, const double& time,
               FunctionType bcFunction)
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
            subRhs.reset(new Vector(curLocalMap));
            subRhs->zero();
            subRhs->subset(*rhs, curLocalMap, offset + suboffset, 0);
            rhss.push_back(subRhs);
            suboffset += curLocalMap.mapSize();
        }
        // apply bcs
        AssemblerType& curAssembler = *it->second;
        (curAssembler.*bcFunction)(rhss, coeff, time);

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
setLawInflow(std::function<double(double)> maxLaw)
{
    typedef std::pair<unsigned int, AssemblerTypePtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;

    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        it->second->setLawInflow(maxLaw);
    }
}

template<class AssemblerType>
void
GlobalAssembler<AssemblerType>::
setLawDtInflow(std::function<double(double)> maxLawDt)
{
    typedef std::pair<unsigned int, AssemblerTypePtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;

    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        it->second->setLawDtInflow(maxLawDt);
    }
}

template<class AssemblerType>
void
GlobalAssembler<AssemblerType>::
setExactSolution(AbstractFunctor* exactSolution)
{
    typedef std::pair<unsigned int, AssemblerTypePtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;

    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        it->second->setExactSolution(exactSolution);
    }
}

template<class AssemblerType>
void
GlobalAssembler<AssemblerType>::
setForcingFunction(Function forcingFunction, Function forcingFunctionDt)
{
    typedef std::pair<unsigned int, AssemblerTypePtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;

    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        it->second->setForcingFunction(forcingFunction, forcingFunctionDt);
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
                                 M_offsets[M_nPrimalBlocks + indices[in]], 0);
            localSolutions.push_back(subSolution);
            in++;
        }
        it->second->exportSolutions(time, localSolutions);
    }
}

template<class AssemblerType>
void
GlobalAssembler<AssemblerType>::
appendNormsToFile(const double& time, VectorPtr solution,
                  std::ofstream& outFile)
{
    typedef std::pair<unsigned int, AssemblerTypePtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;
    typedef std::shared_ptr<LifeV::MapEpetra>            MapEpetraPtr;
    typedef std::vector<MapEpetraPtr>                    MapVector;

    std::vector<double> norms;
    unsigned int countDomains = 0;
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
        std::vector<double> localNorms = it->second->computeNorms(localSolutions);
        if (countDomains == 0)
        {
            for (std::vector<double>::iterator itNorm = localNorms.begin();
                 itNorm != localNorms.end(); itNorm++)
            {
                norms.push_back(*itNorm);
            }
        }
        else
        {
            unsigned int count = 0;
            for (std::vector<double>::iterator itNorm = localNorms.begin();
                 itNorm != localNorms.end(); itNorm++)
            {
                double curNorm = *itNorm;
                double newNorm = curNorm * curNorm +
                                 norms[count] * norms[count];
                norms[count] = std::sqrt(newNorm);
                count++;
            }
        }
        countDomains++;
    }

    std::string newLine = std::to_string(time);
    for (std::vector<double>::iterator itNorm = norms.begin();
         itNorm != norms.end(); itNorm++)
    {
        newLine += ",";
        std::ostringstream streamOb;
        streamOb << *itNorm;
        newLine += streamOb.str();    }
    newLine += "\n";
    outFile << newLine << std::flush;
}

template<class AssemblerType>
void
GlobalAssembler<AssemblerType>::
appendErrorsToFile(const double& time, VectorPtr solution,
                   std::ofstream& outFile)
{
    typedef std::pair<unsigned int, AssemblerTypePtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;
    typedef std::shared_ptr<LifeV::MapEpetra>            MapEpetraPtr;
    typedef std::vector<MapEpetraPtr>                    MapVector;

    std::vector<double> errors;
    unsigned int countDomains = 0;
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
        std::vector<double> localErrors = it->second->computeErrors(localSolutions,
                                                                    time);
        if (countDomains == 0)
        {
            for (std::vector<double>::iterator itError = localErrors.begin();
                 itError != localErrors.end(); itError++)
            {
                errors.push_back(*itError);
            }
        }
        else
        {
            unsigned int count = 0;
            for (std::vector<double>::iterator itErrors = localErrors.begin();
                 itErrors != localErrors.end(); itErrors++)
            {
                double curError = *itErrors;
                double newError = curError * curError +
                                  errors[count] * errors[count];
                errors[count] = std::sqrt(newError);
                count++;
            }
        }
        countDomains++;
    }

    std::string newLine = std::to_string(time);
    for (std::vector<double>::iterator itErrors = errors.begin();
         itErrors != errors.end(); itErrors++)
    {
        newLine += ",";
        std::ostringstream streamOb;
        streamOb << *itErrors;
        newLine += streamOb.str();
    }
    newLine += "\n";
    outFile << newLine << std::flush;
}

template<class AssemblerType>
void
GlobalAssembler<AssemblerType>::
setTimeIntegrationOrder(unsigned int order)
{
    typedef std::pair<unsigned int, AssemblerTypePtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;

    unsigned int offset = 0;
    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        it->second->setTimeIntegrationOrder(order);
    }
}

template<class AssemblerType>
void
GlobalAssembler<AssemblerType>::
checkResidual(VectorPtr solution, VectorPtr prevSolution, double dt)
{
    typedef std::pair<unsigned int, AssemblerTypePtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;
    typedef std::shared_ptr<LifeV::MapEpetra>            MapEpetraPtr;
    typedef std::vector<MapEpetraPtr>                    MapVector;

    unsigned int offsetPrimal = 0;
    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        std::vector<VectorPtr> localSolutions;
        std::vector<VectorPtr> localPrevSolutions;
        // first handle the primal solutions
        MapVector maps = it->second->getPrimalMapVector();
        for (MapVector::iterator itmap = maps.begin();
             itmap != maps.end(); itmap++)
        {
            LifeV::MapEpetra& curLocalMap = **itmap;
            VectorPtr subSolution;
            subSolution.reset(new Vector(curLocalMap));
            subSolution->zero();
            subSolution->subset(*solution, curLocalMap, offsetPrimal, 0);
            localSolutions.push_back(subSolution);
            subSolution.reset(new Vector(curLocalMap));
            subSolution->zero();
            subSolution->subset(*prevSolution, curLocalMap, offsetPrimal, 0);
            localPrevSolutions.push_back(subSolution);
            offsetPrimal += curLocalMap.mapSize();
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
                                 M_offsets[M_nPrimalBlocks + indices[in]], 0);
            localSolutions.push_back(subSolution);
            subSolution.reset(new Vector(curLocalMap));
            subSolution->zero();
            subSolution->subset(*prevSolution, curLocalMap,
                                 M_offsets[M_nPrimalBlocks + indices[in]], 0);
            localPrevSolutions.push_back(subSolution);
            in++;
        }

        it->second->checkResidual(localSolutions,localPrevSolutions,dt);
    }
}

template<class AssemblerType>
void
GlobalAssembler<AssemblerType>::
setTimestep(double dt)
{
    typedef std::pair<unsigned int, AssemblerTypePtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;

    unsigned int offset = 0;
    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        it->second->setTimestep(dt);
    }
}

template<class AssemblerType>
void
GlobalAssembler<AssemblerType>::
postProcess()
{
    typedef std::pair<unsigned int, AssemblerTypePtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;

    unsigned int offset = 0;
    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        it->second->postProcess();
    }
}

template<class AssemblerType>
void
GlobalAssembler<AssemblerType>::
printMeshSize(std::string filename)
{
    typedef std::pair<unsigned int, AssemblerTypePtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;

    std::ofstream outFile;

    outFile.open(filename);
    outFile << "h\n";
    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        outFile << it->second->getMeshSize() << "\n";
    }
    outFile.close();
}

}  // namespace RedMA
