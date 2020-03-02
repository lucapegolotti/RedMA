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
        resize(newSnapshot.nRows(), newSnapshot.nCols());

        for (unsigned int i = 0; i < M_nRows; i++)
        {
            for (unsigned int j = 0; j < M_nCols; j++)
            {
                M_mdeims(i,j)->initialize(newSnapshot.block(i,j));
                M_structures(i,j) = M_mdeims(i,j)->getMDEIMStructure();
            }
        }
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
resize(unsigned int rows, unsigned int cols)
{
    M_nRows = rows;
    M_nCols = cols;
    M_mdeims.resize(M_nRows, M_nCols);
    M_structures.resize(M_nRows, M_nCols);
    setupMDEIMs();
}

void
BlockMDEIM::
setRBBases(SHP(RBBases) bases)
{
    M_bases = bases;
}

void
BlockMDEIM::
setupMDEIMs()
{

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            M_mdeims(i,j).reset(new MDEIM());
            M_mdeims(i,j)->setComm(M_comm);
            M_mdeims(i,j)->setDataContainer(M_data);
            M_mdeims(i,j)->setFESpace(M_assembler->getFEspace(i));
            M_mdeims(i,j)->setRangeMap(M_assembler->getFEspace(i)->mapPtr());
            M_mdeims(i,j)->setDomainMap(M_assembler->getFEspace(j)->mapPtr());
        }
    }
}

void
BlockMDEIM::
performMDEIM(std::string dir)
{
    using namespace boost::filesystem;

    std::string outdir = dir + "/blockmdeim" + std::to_string(M_matIndex);
    create_directory(outdir);

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            std::string curdir = outdir + "/mdeim_" + std::to_string(i) +
                                 "_" + std::to_string(j);

            create_directory(curdir);

            M_mdeims(i,j)->performMDEIM(curdir);
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
    BlockMatrix<FEMATRIX> completeMat = M_assembler->assembleMatrix(M_matIndex);
    BlockMatrix<FEMATRIX> approxMat = assembleMatrix();

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            FEMATRIX curMat(approxMat.block(i,j));
            if (curMat.data())
            {
                std::cout << "=======CHECKING ONLINE MDEIM========" << std::endl << std::flush;
                std::cout << "NORM apprMatrix = " << curMat.data()->normFrobenius() << std::endl << std::flush;

                curMat -= completeMat.block(i,j);

                std::cout << "NORM actualMatrix = " << completeMat.block(i,j).data()->normFrobenius() << std::endl << std::flush;
                std::cout << "REL NORM DIFFERENCE = " << curMat.data()->normFrobenius() / completeMat.block(i,j).data()->normFrobenius() << std::endl << std::flush;
            }
        }
    }
}

BlockMatrix<MatrixEp>
BlockMDEIM::
assembleMatrix()
{
    BlockMatrix<FEMATRIX> retMat;
    retMat.resize(M_nRows, M_nCols);

    BlockMatrix<FEMATRIX> reducedMat = M_assembler->assembleMatrix(M_matIndex, &M_structures);

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            retMat.block(i,j).softCopy(M_mdeims(i,j)->assembleMatrix(reducedMat.block(i,j)));
        }
    }
    return retMat;
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
            M_mdeims(i,j)->projectMDEIM(M_bases->getBasis(i), M_bases->getBasis(j));
        }
    }
}

void
BlockMDEIM::
loadMDEIM(std::string dir)
{
    using namespace boost::filesystem;

    unsigned countrows = 0;
    while (exists(dir + "/mdeim_" + std::to_string(countrows) + "_0"))
        countrows++;

    unsigned countcols = 0;
    while (exists(dir + "/mdeim_0_" + std::to_string(countcols)))
        countcols++;

    resize(countrows, countcols);

    for (unsigned int i = 0; i < M_nRows; i++)
    {
        for (unsigned int j = 0; j < M_nCols; j++)
        {
            std::string curdir = dir + "/mdeim_" + std::to_string(i) + "_" + std::to_string(j);
            M_mdeims(i,j)->loadMDEIM(curdir);
            // we need to redo this because now the structure have changed
            M_structures(i,j) = M_mdeims(i,j)->getMDEIMStructure();
        }
    }
}

}  // namespace RedMA
