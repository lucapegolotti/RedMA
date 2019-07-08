#include <GlobalBlockMatrix.hpp>

namespace RedMA
{

GlobalBlockMatrix::
GlobalBlockMatrix()
{
}

GlobalBlockMatrix::
GlobalBlockMatrix(unsigned int numRows, unsigned int numCols)
{
    resize(numRows, numCols);
}

void
GlobalBlockMatrix::
resize(unsigned int numRows, unsigned int numCols)
{
    std::cout << "resize" << std::endl << std::flush;
    M_rows = numRows;
    M_cols = numCols;
    M_gridEpetra.resize(numRows, numCols);
    for (unsigned int i = 0; i < numRows; i++)
    {
        for (unsigned int j = 0; j < numCols; j++)
        {
            M_gridEpetra(i,j) = nullptr;
        }
    }
}

GlobalBlockMatrix::
GlobalBlockMatrix(const GlobalBlockMatrix& other)
{
    std::cout << "copy" << std::endl << std::flush;
    resize(other.M_rows, other.M_cols);

    for (unsigned int i = 0; i < M_rows; i++)
    {
        for (unsigned int j = 0; j < M_cols; j++)
        {
            if (other.M_gridEpetra(i,j))
            {
                MatrixEpetraPtr newMatrix(new
                                        MatrixEpetra(*other.M_gridEpetra(i,j)));
                M_gridEpetra(i,j) = newMatrix;
            }
        }
    }
    M_globalRangeMap = other.M_globalRangeMap;
    M_globalDomainMap = other.M_globalDomainMap;
    M_domainMaps = other.M_domainMaps;
    M_rangeMaps = other.M_rangeMaps;
}

GlobalBlockMatrix::MatrixEpetraPtr&
GlobalBlockMatrix::
block(unsigned int row, unsigned int col)
{
    std::cout << "block" << std::endl << std::flush;
    return M_gridEpetra(row, col);
}

void
GlobalBlockMatrix::
add(const GlobalBlockMatrix& other)
{
    std::cout << "add" << std::endl << std::flush;

    if (M_rows != other.M_rows && M_cols != other.M_cols)
    {
        throw Exception("Dimensions must be consistent in sum of matrices!");
    }

    for (unsigned int i = 0; i < M_rows; i++)
    {
        for (unsigned int j = 0; j < M_cols; j++)
        {
            if (other.M_gridEpetra(i,j) != nullptr && M_gridEpetra(i,j) != nullptr)
                *M_gridEpetra(i,j) += *other.M_gridEpetra(i,j);
            else if (M_gridEpetra(i,j) == nullptr && M_gridEpetra(i,j) != nullptr)
                M_gridEpetra(i,j).reset(new MatrixEpetra(*other.M_gridEpetra(i,j)));
        }
    }
}

// void
// GlobalBlockMatrix::
// multiply(const GlobalBlockMatrix& other, GlobalBlockMatrix& result)
// {
//     if (M_cols != other.M_rows)
//     {
//         throw Exception("Inner dimensions must be the same in multiplication!");
//     }
//
//     result.resize(M_rows, other.M_cols);
//
//     for (unsigned int i = 0; i < M_rows; i++)
//     {
//         for (unsigned int j = 0; j < other.M_cols; j++)
//         {
//             MatrixEpetraPtr newMatrix = nullptr;
//             for (unsigned int k = 0; k < M_cols; k++)
//             {
//                 if (M_gridEpetra(i,j) != nullptr && other.M_gridEpetra(i,j) != nullptr)
//                 {
//                     if (newMatrix == nullptr)
//                     {
//                         newMatrix.reset(new MatrixEpetra(*M_gridEpetra(i,k)));
//                         *newMatrix *= *other.M_gridEpetra(k,j);
//                     }
//                     else
//                         *newMatrix += (*M_gridEpetra(i,k)) * (*other.M_gridEpetra(k,j));
//                 }
//             }
//             result(i,j) = newMatrix;
//         }
//     }
// }

GlobalBlockMatrix&
GlobalBlockMatrix::
operator*=(const double& coeff)
{
    std::cout << "+=" << std::endl << std::flush;

    for (unsigned int i = 0; i < M_rows; i++)
    {
        for (unsigned int j = 0; j < M_cols; j++)
        {
            if (M_gridEpetra(i,j) != nullptr)
                *M_gridEpetra(i,j) *= coeff;
        }
    }
    return *this;
}

GlobalBlockMatrix::Grid&
GlobalBlockMatrix::
getGrid()
{
    std::cout << "getGrid" << std::endl << std::flush;

    M_grid.resize(M_rows, M_cols);

    for (unsigned int i = 0; i < M_rows; i++)
    {
        for (unsigned int j = 0; j < M_cols; j++)
        {
            if (M_gridEpetra(i,j) == nullptr)
                M_grid(i,j) = nullptr;
            else
                M_grid(i,j) = M_gridEpetra(i,j)->matrixPtr();
        }
    }
    return M_grid;
}

void
GlobalBlockMatrix::
copyBlock(unsigned int rows, unsigned int cols, MatrixEpetraPtr matrix)
{
    std::cout << "copyBlock" << std::endl << std::flush;
    std::cout << "rows = " << rows << "/" << M_rows << std::endl;
    std::cout << "cols = " << cols << "/" << M_cols << std::endl;
    if (matrix)
        M_gridEpetra(rows,cols).reset(new MatrixEpetra(*matrix));
}

// Attention: we assume that the global matrix is square (i.e. dimensions of
// input and output vectors are the same)
GlobalBlockMatrix::VectorEpetra
GlobalBlockMatrix::
operator*(const VectorEpetra& vector)
{
    std::cout << "operator*" << std::endl << std::flush;

    VectorEpetra newVector(*M_globalDomainMap);
    newVector.zero();

    std::vector<VectorEpetraPtr> rangeVectors;
    std::vector<unsigned int> offsets;

    offsets.push_back(0);
    for (unsigned int i = 0; i < M_cols; i++)
    {
        std::cout << "i = " << i << std::endl;
        VectorEpetraPtr rangeVector(new VectorEpetra(*M_rangeMaps[i]));
        rangeVector->subset(vector, *M_rangeMaps[i], offsets[i], 0);
        rangeVectors.push_back(rangeVector);
        offsets.push_back(offsets[i] + M_rangeMaps[i]->mapSize());
    }

    unsigned int offset = 0;
    for (unsigned int i = 0; i < M_rows; i++)
    {
        std::cout << "i = " << i << std::endl;
        VectorEpetraPtr subVector(new VectorEpetra(*M_domainMaps[i]));
        subVector->zero();
        for (unsigned int j = 0; j < M_cols; j++)
        {
            std::cout << "j = " << i << std::endl;
            if (M_gridEpetra(i,j) != nullptr)
            {
                *subVector += (*M_gridEpetra(i,j)) * (*rangeVectors[j]);
            }
        }
        newVector.subset(*subVector, *M_domainMaps[i], 0, offset);
        offset += M_domainMaps[i]->mapSize();
    }

    return newVector;
}

void
GlobalBlockMatrix::
setMaps(std::vector<MapPtr> rangeMaps, std::vector<MapPtr> domainMaps)
{
    std::cout << "createMaps" << std::endl << std::flush;

    M_rangeMaps = rangeMaps;
    M_domainMaps = domainMaps;

    M_globalRangeMap.reset(new Map());
    for (std::vector<MapPtr>::iterator it = M_rangeMaps.begin();
         it != M_rangeMaps.end(); it++)
    {
        if (*it != nullptr)
            *M_globalRangeMap += *(*it);
        else
            throw Exception("Not all range maps were filled!");
    }

    M_globalDomainMap.reset(new Map());
    for (std::vector<MapPtr>::iterator it = M_domainMaps.begin();
         it != M_domainMaps.end(); it++)
    {
        if (*it != nullptr)
            *M_globalDomainMap += *(*it);
        else
            throw Exception("Not all domain maps were filled!");
    }

}


} // namespace RedMA
