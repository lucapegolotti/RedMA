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
addSnapshot(BlockMatrix newSnapshot)
{
    // for (unsigned int i = 0; i < M_nRows; i++)
    // {
    //     for (unsigned int j = 0; j < M_nCols; j++)
    //     {
    //         M_mdeims(i,j)->addSnapshot(newSnapshot.block(i,j));
    //     }
    // }
}

void
BlockMDEIM::
setFESpace(shp<FESPACE> fespace, const unsigned int& index)
{
    // M_fespaces[index] = fespace;
    // for (unsigned int j = 0; j < M_nCols; j++)
    // {
    //     if (M_mdeims(index,j))
    //     {
    //         M_mdeims(index,j)->setFESpace(fespace);
    //     }
    // }
}

void
BlockMDEIM::
setRangeMap(shp<MAPEPETRA> map, const unsigned int& index)
{
    // M_rangeMaps[index] = map;
    //
    // for (unsigned int j = 0; j < M_nCols; j++)
    // {
    //     if (M_mdeims(index,j))
    //         M_mdeims(index,j)->setRangeMap(map);
    // }
}

void
BlockMDEIM::
setDomainMap(shp<MAPEPETRA> map, const unsigned int& index)
{
    // M_domainMaps[index] = map;
    //
    // for (unsigned int i = 0; i < M_nRows; i++)
    // {
    //     if (M_mdeims(i,index))
    //         M_mdeims(i,index)->setDomainMap(map);
    // }
}

void
BlockMDEIM::
resize(unsigned int rows, unsigned int cols)
{
    // M_nRows = rows;
    // M_nCols = cols;
    // M_mdeims.resize(M_nRows, M_nCols);
    // M_fespaces.resize(M_nRows);
    // M_rangeMaps.resize(M_nRows);
    // M_domainMaps.resize(M_nCols);
    // M_structures.resize(M_nRows, M_nCols);
    //
    // for (unsigned int i = 0; i < M_nRows; i++)
    // {
    //     for (unsigned int j = 0; j < M_nCols; j++)
    //     {
    //         M_mdeims(i,j).reset(new MDEIM());
    //     }
    // }
}

void
BlockMDEIM::
setRBBases(shp<RBBases> bases)
{
    M_bases = bases;
}

void
BlockMDEIM::
performMDEIM(std::string dir)
{
    // // using namespace boost::filesystem;
    //
    // std::string outdir = dir + "/blockmdeim" + std::to_string(M_matIndex);
    // create_directory(outdir);
    //
    // for (unsigned int i = 0; i < M_nRows; i++)
    // {
    //     for (unsigned int j = 0; j < M_nCols; j++)
    //     {
    //         std::string curdir = outdir + "/mdeim_" + std::to_string(i) +
    //                              "_" + std::to_string(j);
    //
    //         create_directory(curdir);
    //
    //         M_mdeims(i,j)->performMDEIM(curdir);
    //     }
    // }
}

void
BlockMDEIM::
prepareOnline()
{
    // for (unsigned int i = 0; i < M_nRows; i++)
    // {
    //     for (unsigned int j = 0; j < M_nCols; j++)
    //     {
    //         M_mdeims(i,j)->prepareOnline();
    //     }
    // }
}

void
BlockMDEIM::
initialize(BlockMatrix reducedMatrix)
{
    // for (unsigned int i = 0; i < M_nRows; i++)
    // {
    //     for (unsigned int j = 0; j < M_nCols; j++)
    //     {
    //         M_mdeims(i,j)->setComm(M_comm);
    //         M_mdeims(i,j)->setDataContainer(M_data);
    //         M_mdeims(i,j)->initialize(reducedMatrix.block(i,j));
    //         M_structures(i,j) = M_mdeims(i,j)->getMDEIMStructure();
    //     }
    // }
}

void
BlockMDEIM::
checkOnline(BlockMatrix completeMat, BlockMatrix reducedMat)
{
    // BlockMatrix approxMat = assembleMatrix(reducedMat);
    //
    // for (unsigned int i = 0; i < M_nRows; i++)
    // {
    //     for (unsigned int j = 0; j < M_nCols; j++)
    //     {
    //         FEMATRIX curMat(approxMat.block(i,j));
    //         if (curMat.data())
    //         {
    //             std::cout << "=======CHECKING ONLINE MDEIM========" << std::endl << std::flush;
    //             std::cout << "NORM apprMatrix = " << curMat.data()->normFrobenius() << std::endl << std::flush;
    //
    //             curMat -= completeMat.block(i,j);
    //
    //             std::cout << "NORM actualMatrix = " << completeMat.block(i,j).data()->normFrobenius() << std::endl << std::flush;
    //             std::cout << "REL NORM DIFFERENCE = " << curMat.data()->normFrobenius() / completeMat.block(i,j).data()->normFrobenius() << std::endl << std::flush;
    //         }
    //     }
    // }
}

BlockMatrix
BlockMDEIM::
assembleMatrix(BlockMatrix reducedMat)
{
    // BlockMatrix retMat;
    // retMat.resize(M_nRows, M_nCols);
    //
    // for (unsigned int i = 0; i < M_nRows; i++)
    // {
    //     for (unsigned int j = 0; j < M_nCols; j++)
    //     {
    //         retMat.block(i,j).shallowCopy(M_mdeims(i,j)->assembleMatrix(reducedMat.block(i,j)));
    //     }
    // }
    // return retMat;
}

BlockMatrix
BlockMDEIM::
assembleProjectedMatrix(BlockMatrix reducedMat)
{
    // BlockMatrix retMat;
    // retMat.resize(M_nRows, M_nCols);
    //
    // for (unsigned int i = 0; i < M_nRows; i++)
    // {
    //     for (unsigned int j = 0; j < M_nCols; j++)
    //     {
    //         retMat.block(i,j).shallowCopy(M_mdeims(i,j)->assembleProjectedMatrix(reducedMat.block(i,j)));
    //     }
    // }
    // return retMat;
}

void
BlockMDEIM::
checkOnSnapshots()
{
    // for (unsigned int i = 0; i < M_nRows; i++)
    // {
    //     for (unsigned int j = 0; j < M_nCols; j++)
    //     {
    //         M_mdeims(i,j)->checkOnSnapshots();
    //     }
    // }
}

void
BlockMDEIM::
dumpMDEIMs(std::string dir, std::string prefix)
{
    // // using namespace boost::filesystem;
    //
    // std::string outdir = dir + "/" + prefix + std::to_string(M_matIndex);
    // create_directory(outdir);
    //
    // for (unsigned int i = 0; i < M_nRows; i++)
    // {
    //     for (unsigned int j = 0; j < M_nCols; j++)
    //     {
    //         std::string curdir = outdir + "/mdeim_" + std::to_string(i) + "_" + std::to_string(j);
    //
    //         create_directory(curdir);
    //         M_mdeims(i,j)->dumpMDEIM(curdir);
    //     }
    // }
}

void
BlockMDEIM::
projectMDEIMs()
{
    // for (unsigned int i = 0; i < M_nRows; i++)
    // {
    //     for (unsigned int j = 0; j < M_nCols; j++)
    //     {
    //         M_mdeims(i,j)->projectMDEIM(M_bases->getEnrichedBasis(i, 0),
    //                                     M_bases->getEnrichedBasis(j, 0));
    //     }
    // }
}

void
BlockMDEIM::
loadMDEIM(std::string dir)
{
    // // using namespace boost::filesystem;
    //
    // for (unsigned int i = 0; i < M_nRows; i++)
    // {
    //     for (unsigned int j = 0; j < M_nCols; j++)
    //     {
    //         std::string curdir = dir + "/mdeim_" + std::to_string(i) + "_" + std::to_string(j);
    //
    //         M_mdeims(i,j)->setComm(M_comm);
    //         M_mdeims(i,j)->setDataContainer(M_data);
    //         M_mdeims(i,j)->loadMDEIM(curdir);
    //         M_structures(i,j) = M_mdeims(i,j)->getMDEIMStructure();
    //     }
    // }
}

}  // namespace RedMA
