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

        newAssembler->addMaps(M_globalMap, M_dimensionsVector);
        M_assemblersVector.push_back(std::make_pair(it->first,newAssembler));
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
getJacobianF()
{
    MatrixPtr jacobian;
    fillGlobalMatrix(jacobian, &AssemblerType::getJacobian);
    return jacobian;
}

template <class AssemblerType>
typename GlobalAssembler<AssemblerType>::VectorPtr
GlobalAssembler<AssemblerType>::
computeF(const double& time, VectorPtr u) const
{
    VectorPtr f(new Vector(*M_globalMap));

    // implementation here
    return f;
}

template <class AssemblerType>
typename GlobalAssembler<AssemblerType>::VectorPtr
GlobalAssembler<AssemblerType>::
computeFder(const double& time, VectorPtr u) const
{
    VectorPtr fder(new Vector(*M_globalMap));

    // implementation here
    return fder;
}

template <class AssemblerType>
void
GlobalAssembler<AssemblerType>::
assembleGlobalMass()
{
    fillGlobalMatrix(M_massMatrix, &AssemblerType::getMassMatrix);
}

template<class AssemblerType>
template<typename FunctionType>
void
GlobalAssembler<AssemblerType>::
fillGlobalMatrix(MatrixPtr& matrixToFill, FunctionType getMatrixMethod)
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
    matrixToFill->globalAssemble();
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
    for (typename AssemblersVector::iterator it = M_assemblersVector.begin();
         it != M_assemblersVector.end(); it++)
    {
        std::vector<VectorPtr> localSolutions;
        MapVector maps = it->second->getMapVector();
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
        it->second->setTimeAndPrevSolution(time, localSolutions);
    }
}

}  // namespace RedMA
