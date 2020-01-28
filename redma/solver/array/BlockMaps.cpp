#include "BlockMap.hpp"

namespace RedMA
{

class BlockMap
{

BlockMap::
BlockMap(const MatrixType& matrix) :
  M_matrix(matrix)
{
    generateMaps();
}

void generateMaps()
{
    unsigned int nRows = M_matrix.nRows();
    unsigned int nCols = M_matrix.nCols();
}

}
