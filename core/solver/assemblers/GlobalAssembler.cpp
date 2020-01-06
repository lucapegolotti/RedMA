#include <GlobalAssembler.hpp>

namespace RedMA
{

GlobalAssembler::
GlobalAssembler(const GetPot& datafile, commPtr_Type comm, bool verbose) :
  AbstractAssembler(datafile,comm,nullptr,verbose),
  M_datafile(datafile),
  M_comm(comm),
  M_verbose(verbose)
{
    M_globalMap.reset(new MapEpetra());
}

void
GlobalAssembler::
buildPrimalStructures(TreeStructure& tree)
{
    typedef std::map<unsigned int, TreeNodePtr>     NodesMap;

    NodesMap nodesMap = tree.getNodesMap();

    for (NodesMap::iterator it = nodesMap.begin(); it != nodesMap.end(); it++)
    {
        AbstractAssemblerPtr newAssembler = AssemblersFactory(M_datafile, M_comm,
                                                              it->second, M_verbose);
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

unsigned int
GlobalAssembler::
numberOfComponents()
{
    throw new Exception("numberOfComponents does not make sense for GlobalAssembler!");
    return 0;
}

void
GlobalAssembler::
setTreeStructure(TreeStructure tree)
{
    M_tree = tree;
}

void
GlobalAssembler::
setup()
{
    typedef std::vector<std::pair<unsigned int, AbstractAssemblerPtr> > AssemblersVector;

    buildPrimalStructures(M_tree);
    buildDualStructures(M_tree);

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

void
GlobalAssembler::
buildDualStructures(TreeStructure& tree)
{
    typedef std::map<unsigned int, TreeNodePtr>                         NodesMap;
    typedef std::vector<TreeNodePtr>                                    NodesVector;
    typedef std::vector<std::pair<unsigned int, AbstractAssemblerPtr> > AssemblersVector;

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
        AbstractAssemblerPtr fatherAssembler;

        fatherAssembler = M_assemblersMap[myID];
        for (NodesVector::iterator itVector = children.begin();
             itVector != children.end(); itVector++)
        {
            if (*itVector)
            {
                unsigned int otherID = (*itVector)->M_ID;
                M_interfaces.push_back(std::make_pair(myID, otherID));

                AbstractAssemblerPtr childAssembler;
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

GlobalAssembler::MapEpetraPtr
GlobalAssembler::
getGlobalMap() const
{
    return M_globalMap;
}

BlockMatrix
GlobalAssembler::
getGlobalMass()
{
    BlockMatrix updatedMassMatrix(M_massMatrix.getNumberRows(),
                                  M_massMatrix.getNumberCols());
    double diagCoefficient = 0.0;
    fillGlobalMatrix(updatedMassMatrix, false, &AbstractAssembler::getUpdateMass,
                     &diagCoefficient);
    updatedMassMatrix.add(M_massMatrix);
    return updatedMassMatrix;
    // return M_massMatrix;
}

GlobalAssembler::VectorPtr
GlobalAssembler::
getInitialCondition()
{
    VectorPtr f(new Vector(*M_globalMap));
    fillGlobalVector(f, &AbstractAssembler::initialCondition);
    return f;
}

BlockMatrix
GlobalAssembler::
getGlobalMassJac()
{
    BlockMatrix updatedMassMatrix(M_massMatrix.getNumberRows(),
                                  M_massMatrix.getNumberCols());
    double diagCoefficient = 0.0;
    fillGlobalMatrix(updatedMassMatrix, false, &AbstractAssembler::getUpdateMassJac,
                     &diagCoefficient);
    updatedMassMatrix.add(M_massMatrix);
    // return updatedMassMatrix;
    return updatedMassMatrix;
}

BlockMatrix
GlobalAssembler::
getGlobalMassJacVelocity()
{
    BlockMatrix updatedMassMatrix(M_massMatrix.getNumberRows(),
                                  M_massMatrix.getNumberCols());
    double diagCoefficient = 0.0;
    fillGlobalMatrix(updatedMassMatrix, false, &AbstractAssembler::getUpdateMassJacVelocity,
                     &diagCoefficient);
    return updatedMassMatrix;
}

BlockMatrix
GlobalAssembler::
getJacobianF(bool addCoupling, double* diagonalCoefficient)
{
    BlockMatrix jacobian;
    fillGlobalMatrix(jacobian, addCoupling, &AbstractAssembler::getJacobian,
                     diagonalCoefficient);
    return jacobian;
}

BlockMatrix
GlobalAssembler::
getJacobianFprec(bool addCoupling, double* diagonalCoefficient)
{
    BlockMatrix jacobian;
    fillGlobalMatrix(jacobian, addCoupling, &AbstractAssembler::getJacobianPrec,
                     diagonalCoefficient);
    return jacobian;
}

GlobalAssembler::VectorPtr
GlobalAssembler::
computeF_()
{
    VectorPtr f(new Vector(*M_globalMap));
    fillGlobalVector(f, &AbstractAssembler::computeF);
    return f;
}

GlobalAssembler::VectorPtr
GlobalAssembler::
computeFder_()
{
    VectorPtr fder(new Vector(*M_globalMap));
    fillGlobalVector(fder, &AbstractAssembler::computeFder);
    return fder;
}

BlockMatrix
GlobalAssembler::
assembleGlobalMass(bool addCoupling, double* diagonalCoefficient)
{
    fillGlobalMatrix(M_massMatrix, addCoupling, &AbstractAssembler::getMassMatrix,
                     diagonalCoefficient);
    return M_massMatrix;
}

void
GlobalAssembler::
setTimeAndPrevSolution(const double& time, VectorPtr solution, bool doAssembly)
{
    typedef std::pair<unsigned int, AbstractAssemblerPtr>    Pair;
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

void
GlobalAssembler::
applyBCsRhsRosenbrock(VectorPtr rhs, VectorPtr utilde,
                      const double& time, const double& dt,
                      const double& alphai, const double& gammai)
{
    typedef std::pair<unsigned int, AbstractAssemblerPtr>    Pair;
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

void
GlobalAssembler::
setLawInflow(std::function<double(double)> maxLaw)
{
    typedef std::pair<unsigned int, AbstractAssemblerPtr>    Pair;
    typedef std::vector<Pair>                                AssemblersVector;

    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        it->second->setLawInflow(maxLaw);
    }
}

void
GlobalAssembler::
setLawDtInflow(std::function<double(double)> maxLawDt)
{
    typedef std::pair<unsigned int, AbstractAssemblerPtr>    Pair;
    typedef std::vector<Pair>                                AssemblersVector;

    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        it->second->setLawDtInflow(maxLawDt);
    }
}

void
GlobalAssembler::
setExactSolution(AbstractFunctor* exactSolution)
{
    typedef std::pair<unsigned int, AbstractAssemblerPtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;

    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        it->second->setExactSolution(exactSolution);
    }
}

void
GlobalAssembler::
setForcingFunction(Function forcingFunction, Function forcingFunctionDt)
{
    typedef std::pair<unsigned int, AbstractAssemblerPtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;

    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        it->second->setForcingFunction(forcingFunction, forcingFunctionDt);
    }
}

void
GlobalAssembler::
exportSolutions(const double& time, VectorPtr solution)
{
    typedef std::pair<unsigned int, AbstractAssemblerPtr>    Pair;
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

void
GlobalAssembler::
appendNormsToFile(const double& time, VectorPtr solution,
                  std::ofstream& outFile)
{
    typedef std::pair<unsigned int, AbstractAssemblerPtr>    Pair;
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

void
GlobalAssembler::
appendErrorsToFile(const double& time, VectorPtr solution,
                   std::ofstream& outFile)
{
    typedef std::pair<unsigned int, AbstractAssemblerPtr>    Pair;
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

void
GlobalAssembler::
setTimeIntegrationOrder(unsigned int order)
{
    typedef std::pair<unsigned int, AbstractAssemblerPtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;

    unsigned int offset = 0;
    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        it->second->setTimeIntegrationOrder(order);
    }
}

void
GlobalAssembler::
setTimestep(double dt)
{
    typedef std::pair<unsigned int, AbstractAssemblerPtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;

    unsigned int offset = 0;
    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        it->second->setTimestep(dt);
    }
}

void
GlobalAssembler::
postProcess()
{
    typedef std::pair<unsigned int, AbstractAssemblerPtr>    Pair;
    typedef std::vector<Pair>                            AssemblersVector;

    unsigned int offset = 0;
    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        it->second->postProcess();
    }
}

void
GlobalAssembler::
printMeshSize(std::string filename)
{
    typedef std::pair<unsigned int, AbstractAssemblerPtr>    Pair;
    typedef std::vector<Pair>                                AssemblersVector;

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
