#include "BlockMDEIM.hpp"

namespace RedMA
{

BlockMDEIM::
BlockMDEIM()
{

}

void
BlockMDEIM::
setDataContainer(const DataContainer& dataContainer)
{
    M_data = dataContainer;
}

void
BlockMDEIM::
setComm(EPETRACOMM comm)
{
    M_comm = comm;
}

void
BlockMDEIM::
setAssembler(SHP(aAssembler<FEVECTOR COMMA FEMATRIX>) assembler)
{
    M_assembler = assembler;
}

// this method works also as a "setup" method for the mdeims
void
BlockMDEIM::
addSnapshot(BlockMatrix<MatrixEp> newSnapshot)
{
    if (M_mdeims.size1() == 0 && M_mdeims.size2() == 0)
    {
        M_nRows = newSnapshot.nRows();
        M_nCols = newSnapshot.nCols();
        M_mdeims.resize(M_nRows, M_nCols);
        M_structures.resize(M_nRows, M_nCols);

        for (unsigned int i = 0; i < M_nRows; i++)
        {
            for (unsigned int j = 0; j < M_nCols; j++)
            {
                M_mdeims(i,j).reset(new MDEIM());
            }
        }
        initializeMDEIMs(newSnapshot);
    }

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            M_mdeims(i,j)->addSnapshot(newSnapshot.block(i,j));
        }
    }
}

void
BlockMDEIM::
setRBBases(std::vector<std::vector<SHP(VECTOREPETRA)>> bases)
{
    M_RBbases = bases;
}

void
BlockMDEIM::
initializeMDEIMs(BlockMatrix<MatrixEp> matrix)
{

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            M_mdeims(i,j)->setComm(M_comm);
            M_mdeims(i,j)->setDataContainer(M_data);
            M_mdeims(i,j)->initialize(matrix.block(i,j));
            M_mdeims(i,j)->setFESpace(M_assembler->getFEspace(i));
            M_mdeims(i,j)->setRangeMap(M_assembler->getFEspace(i)->mapPtr());
            M_mdeims(i,j)->setDomainMap(M_assembler->getFEspace(j)->mapPtr());
            M_structures(i,j) = M_mdeims(i,j)->getMDEIMStructure();
        }
    }
}

void
BlockMDEIM::
performMDEIM()
{
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            M_mdeims(i,j)->performMDEIM();
        }
    }
}

void
BlockMDEIM::
prepareOnline()
{
    BlockMatrix<FEMATRIX> mat = M_assembler->assembleMatrix(M_matIndex, &M_structures);

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            M_mdeims(i,j)->prepareOnline(mat.block(i,j));
        }
    }
}

void
BlockMDEIM::
checkOnline()
{
    BlockMatrix<FEMATRIX> reducedMat = M_assembler->assembleMatrix(M_matIndex, &M_structures);
    BlockMatrix<FEMATRIX> completeMat = M_assembler->assembleMatrix(M_matIndex);

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            M_mdeims(i,j)->checkOnline(reducedMat.block(i,j), completeMat.block(i,j));
        }
    }
}

void
BlockMDEIM::
checkOnSnapshots()
{
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            M_mdeims(i,j)->checkOnSnapshots();
        }
    }
}

void
BlockMDEIM::
dumpMDEIMs(std::string dir)
{
    using namespace boost::filesystem;

    std::string outdir = dir + "/blockmdeim" + std::to_string(M_matIndex);
    create_directory(outdir);

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            std::string curdir = outdir + "/mdeim_" + std::to_string(i) + "_" + std::to_string(j);

            create_directory(curdir);
            M_mdeims(i,j)->dumpMDEIM(curdir);
        }
    }
}

void
BlockMDEIM::
projectMDEIMs()
{
    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            M_mdeims(i,j)->projectMDEIM(M_RBbases[i], M_RBbases[j]);
        }
    }
}

}  // namespace RedMA
