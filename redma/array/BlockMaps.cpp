#include "BlockMaps.hpp"

namespace RedMA
{

// template <>
// void
// BlockMaps<MatrixEp>::
// generateMaps()
// {
//     unsigned int nrows = M_matrix.nRows();
//     unsigned int ncols = M_matrix.nCols();
//
//     M_rangeMaps.resize(nrows);
//     M_domainMaps.resize(ncols);
//     M_rangeMapsEpetra.resize(nrows);
//     M_domainMapsEpetra.resize(ncols);
//     M_nInnerRows.push_back(nrows);
//     M_nInnerCols.push_back(ncols);
//
//     for (unsigned int i = 0; i < nrows; i++)
//     {
//         for (unsigned int j = 0; j < ncols; j++)
//         {
//             M_rangeMaps[i] = M_matrix.block(i,j).data()->rangeMap().map(LifeV::Unique);
//             M_domainMaps[j] = M_matrix.block(i,j).data()->domainMap().map(LifeV::Unique);
//             M_rangeMapsEpetra[i].reset(new MAPEPETRA(M_matrix.block(i,j).data()->rangeMap()));
//             M_domainMapsEpetra[j].reset(new MAPEPETRA(M_matrix.block(i,j).data()->domainMap()));
//         }
//     }
// }
//
// template <>
// void
// BlockMaps<BlockMatrix<MatrixEp>>::
// generateMaps()
// {
//     unsigned int nrows = M_matrix.nRows();
//     unsigned int ncols = M_matrix.nCols();
//
//     // first we deal with ranges
//     for (unsigned int i = 0; i < nrows; i++)
//     {
//         for (unsigned int j = 0; j < ncols; j++)
//         {
//             BlockMaps<MatrixEp> localMaps(M_matrix.block(i,j));
//             auto ranges = localMaps.getRangeMaps();
//
//             for (auto map : ranges)
//                 M_rangeMaps.push_back(map);
//
//             auto rangesEpetra = localMaps.getRangeMapsEpetra();
//             for (auto map : rangesEpetra)
//                 M_rangeMapsEpetra.push_back(map);
//
//             M_nInnerRows.push_back(localMaps.getInnerRows()[0]);
//             break;
//         }
//     }
//
//     for (unsigned int j = 0; j < ncols; j++)
//     {
//         for (unsigned int i = 0; i < nrows; i++)
//         {
//             BlockMaps<MatrixEp > localMaps(M_matrix.block(i,j));
//             auto ranges = localMaps.getDomainMaps();
//
//             for (auto map : ranges)
//                 M_domainMaps.push_back(map);
//
//             auto domainEpetra = localMaps.getDomainMapsEpetra();
//             for (auto map : domainEpetra)
//                 M_domainMapsEpetra.push_back(map);
//
//             M_nInnerCols.push_back(localMaps.getInnerCols()[0]);
//             break;
//         }
//     }
// }
//
// template <>
// BlockMatrix<MatrixEp>
// collapseBlocks(const BlockMatrix<BlockMatrix<MatrixEp>>& matrix,
//                const BlockMaps<BlockMatrix<MatrixEp>>& maps)
// {
//     if (!matrix.isFinalized())
//     {
//         throw new Exception("Matrix must be finalized in order to collapse blocks!");
//     }
//
//     auto rangeMaps = maps.getRangeMaps();
//     auto domainMaps = maps.getRangeMaps();
//
//     unsigned int nrows = rangeMaps.size();
//     unsigned int ncols = domainMaps.size();
//
//     std::vector<unsigned int> ninnerrows = maps.getInnerRows();
//     std::vector<unsigned int> ninnercols = maps.getInnerCols();
//
//     BlockMatrix<MatrixEp> retMatrix;
//     retMatrix.resize(nrows, ncols);
//
//     unsigned int offsetrow = 0;
//     for (unsigned int iouter = 0; iouter < matrix.nRows(); iouter++)
//     {
//         unsigned int offsetcol = 0;
//         for (unsigned int jouter = 0; jouter < matrix.nCols(); jouter++)
//         {
//             auto curmatrix = matrix.block(iouter,jouter);
//
//             for (unsigned int iinner = 0; iinner < ninnerrows[iouter]; iinner++)
//             {
//                 for (unsigned int jinner = 0; jinner < ninnercols[jouter]; jinner++)
//                 {
//                     retMatrix.block(offsetrow+iinner,offsetcol+jinner).
//                               softCopy(curmatrix.block(iinner,jinner));
//                 }
//             }
//             offsetcol += ninnercols[jouter];
//         }
//         offsetrow += ninnerrows[iouter];
//     }
//
//     retMatrix.finalize();
//
//     return retMatrix;
// }
//
// template <>
// MatrixEp
// collapseBlocks(const BlockMatrix<MatrixEp>& matrix,
//                const BlockMaps<MatrixEp>& maps)
// {
//     using namespace LifeV::MatrixEpetraStructuredUtility;
//
//     MatrixEp retMatrix;
//     auto rangeMaps = maps.getRangeMapsEpetra();
//     auto domainMaps = maps.getRangeMapsEpetra();
//
//     retMatrix.data().reset(new MATRIXEPETRA(*maps.getMonolithicRangeMapEpetra()));
//
//     std::vector<unsigned int> dimensionsRows;
//     std::vector<unsigned int> dimensionsCols;
//
//     for (auto map : rangeMaps)
//         dimensionsRows.push_back(map->mapSize());
//
//     for (auto map : domainMaps)
//         dimensionsCols.push_back(map->mapSize());
//
//     LifeV::MatrixBlockStructure structure;
//     structure.setBlockStructure(dimensionsRows,
//                                 dimensionsCols);
//
//     for (unsigned int i = 0; i < matrix.nRows(); i++)
//     {
//         for (unsigned int j = 0; j < matrix.nCols(); j++)
//         {
//             SHP(LifeV::MatrixEpetraStructuredView<double>) globalView;
//             globalView = createBlockView(retMatrix.data(), structure, i, j);
//
//             LifeV::MatrixBlockStructure blockStructure;
//             std::vector<unsigned int> rows(1), cols(1);
//             rows[0] = dimensionsRows[i];
//             cols[0] = dimensionsCols[j];
//             blockStructure.setBlockStructure(rows, cols);
//
//             SHP(LifeV::MatrixEpetraStructuredView<double>) blockLocalView;
//             blockLocalView = createBlockView(matrix.block(i,j).data(),
//                                              blockStructure, 0, 0);
//
//             copyBlock(blockLocalView, globalView);
//         }
//     }
//
//     retMatrix.data()->globalAssemble(maps.getMonolithicDomainMapEpetra(),
//                                      maps.getMonolithicRangeMapEpetra());
//
//     return retMatrix;
// }
//
// SHP(VECTOREPETRA)
// getEpetraVector(const BlockVector<BlockVector<VectorEp>>& vector,
//                 const BlockMaps<BlockMatrix<MatrixEp>>& maps)
// {
//
//     SHP(VECTOREPETRA) retVec(new VECTOREPETRA(maps.getMonolithicRangeMapEpetra(),
//                                               LifeV::Unique));
//     auto rangeMaps = maps.getRangeMapsEpetra();
//     auto rows = maps.getInnerRows();
//
//     unsigned int count = 0;
//     unsigned int offset = 0;
//     // this part is tricky because we don't know a priori if all the
//     // blocks in the input vector are filled
//     for (unsigned int iouter = 0; iouter < vector.nRows(); iouter++)
//     {
//         auto curblock = vector.block(iouter);
//         unsigned int currows = vector.block(iouter).nRows();
//
//         if (curblock.nRows() == 0)
//         {
//             for (unsigned int iinner = 0; iinner < rows[iouter]; iinner++)
//             {
//                 offset += rangeMaps[count]->mapSize();
//                 count++;
//             }
//         }
//         else
//         {
//             for (unsigned int iinner = 0; iinner < curblock.nRows(); iinner++)
//             {
//                 if (curblock.block(iinner).data())
//                 {
//                     retVec->subset(*curblock.block(iinner).data(),
//                                    *rangeMaps[count], 0, offset);
//                 }
//
//                 offset += rangeMaps[count]->mapSize();
//                 count++;
//             }
//
//         }
//
//     }
//
//     return retVec;
// }
//
// BlockVector<BlockVector<VectorEp>>
// getBlockVector(const SHP(VECTOREPETRA)& vector,
//                const BlockMaps<BlockMatrix<MatrixEp>>& maps)
// {
//     BlockVector<BlockVector<VectorEp>> retVec;
//     std::vector<unsigned int> rows = maps.getInnerRows();
//     auto rangeMaps = maps.getRangeMapsEpetra();
//
//     retVec.resize(rows.size());
//
//     unsigned int count = 0;
//     unsigned int offset = 0;
//     for (unsigned int iouter = 0; iouter < rows.size(); iouter++)
//     {
//         retVec.block(iouter).resize(rows[iouter]);
//         for (unsigned int iinner = 0; iinner < rows[iouter]; iinner++)
//         {
//             SHP(VECTOREPETRA)& cvHandler = retVec.block(iouter).block(iinner).data();
//             cvHandler.reset(new VECTOREPETRA(rangeMaps[count], LifeV::Unique));
//             cvHandler->subset(*vector, *rangeMaps[count], offset, 0);
//             offset += rangeMaps[count]->mapSize();
//             count++;
//         }
//     }
//
//     return retVec;
// }

}
