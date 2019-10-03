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

unsigned int
GlobalBlockMatrix::
getNumberRows()
{
    return M_rows;
}

unsigned int
GlobalBlockMatrix::
getNumberCols()
{
    return M_cols;
}

GlobalBlockMatrix::
GlobalBlockMatrix(const GlobalBlockMatrix& other)
{
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
    return M_gridEpetra(row, col);
}

GlobalBlockMatrix::MatrixEpetraPtr
GlobalBlockMatrix::
block(unsigned int row, unsigned int col) const
{
    return M_gridEpetra(row, col);
}

void
GlobalBlockMatrix::
add(const GlobalBlockMatrix& other)
{
    if (M_rows != other.M_rows && M_cols != other.M_cols)
    {
        throw Exception("Dimensions must be consistent in sum of matrices!");
    }

    for (unsigned int i = 0; i < M_rows; i++)
    {
        for (unsigned int j = 0; j < M_cols; j++)
        {
            if (M_gridEpetra(i,j) != nullptr && other.M_gridEpetra(i,j) != nullptr)
            {
                // we do this because if the matrix is zero then we must
                // use the sparsity pattern of the other matrix
                if (M_gridEpetra(i,j)->norm1() > 0 &&
                    other.M_gridEpetra(i,j)->norm1() > 0)
                    *M_gridEpetra(i,j) += *other.M_gridEpetra(i,j);
                else if (M_gridEpetra(i,j)->norm1() == 0)
                    M_gridEpetra(i,j).reset(new MatrixEpetra(*other.M_gridEpetra(i,j)));
            }
            else if (M_gridEpetra(i,j) == nullptr && other.M_gridEpetra(i,j) != nullptr)
            {
                M_gridEpetra(i,j).reset(new MatrixEpetra(*other.M_gridEpetra(i,j)));
            }
        }
    }
    M_globalRangeMap = other.M_globalRangeMap;
    M_globalDomainMap = other.M_globalDomainMap;
    M_domainMaps = other.M_domainMaps;
    M_rangeMaps = other.M_rangeMaps;
}

GlobalBlockMatrix&
GlobalBlockMatrix::
operator*=(const double& coeff)
{
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

GlobalBlockMatrix::Grid
GlobalBlockMatrix::
getGrid()
{
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
    if (matrix)
    {
        M_gridEpetra(rows,cols).reset(new MatrixEpetra(*matrix));
    }

}

// Attention: we assume that the global matrix is square (i.e. dimensions of
// input and output vectors are the same)
GlobalBlockMatrix::VectorEpetra
GlobalBlockMatrix::
operator*(const VectorEpetra& vector)
{
    VectorEpetra newVector(*M_globalDomainMap);
    newVector.zero();

    std::vector<VectorEpetraPtr> rangeVectors;
    std::vector<unsigned int> offsets;

    offsets.push_back(0);
    for (unsigned int i = 0; i < M_cols; i++)
    {
        VectorEpetraPtr rangeVector(new VectorEpetra(*M_rangeMaps[i]));
        rangeVector->subset(vector, *M_rangeMaps[i], offsets[i], 0);
        rangeVectors.push_back(rangeVector);
        offsets.push_back(offsets[i] + M_rangeMaps[i]->mapSize());
    }

    unsigned int offset = 0;
    for (unsigned int i = 0; i < M_rows; i++)
    {
        VectorEpetraPtr subVector(new VectorEpetra(*M_domainMaps[i]));
        subVector->zero();
        for (unsigned int j = 0; j < M_cols; j++)
        {
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

GlobalBlockMatrix::MapPtr
GlobalBlockMatrix::
rangeMap(unsigned int row, unsigned int col) const
{
    if (M_gridEpetra(row,col) != nullptr)
        return std::make_shared<Map>(M_gridEpetra(row,col)->rangeMap());
    return nullptr;
}

GlobalBlockMatrix::MapPtr
GlobalBlockMatrix::
domainMap(unsigned int row, unsigned int col) const
{
    if (M_gridEpetra(row,col) != nullptr)
        return std::make_shared<Map>(M_gridEpetra(row,col)->domainMap());
    return nullptr;
}

void
GlobalBlockMatrix::
printPattern()
{
    for (unsigned int i = 0; i < M_rows; i++)
    {
        for (unsigned int j = 0; j < M_cols; j++)
        {
            if (M_gridEpetra(i,j) != nullptr &&
                M_gridEpetra(i,j)->norm1() > 0)
                std::cout << "x\t";
            else
                std::cout << "o\t";
        }
        std::cout << std::endl;
    }
}

void
GlobalBlockMatrix::
spy()
{
    for (unsigned int i = 0; i < M_rows; i++)
    {
        for (unsigned int j = 0; j < M_cols; j++)
        {
            if (M_gridEpetra(i,j) != nullptr &&
                M_gridEpetra(i,j)->norm1() > 0)
                M_gridEpetra(i,j)->spy("block" + std::to_string(i) + "_"
                                               + std::to_string(j));
        }
    }
}

void
GlobalBlockMatrix::
singleNorms1()
{
    for (unsigned int i = 0; i < M_rows; i++)
    {
        for (unsigned int j = 0; j < M_cols; j++)
        {
            if (M_gridEpetra(i,j) != nullptr &&
                M_gridEpetra(i,j)->norm1() > 0)
            {
                std::cout << "Block(" << std::to_string(i) << ","
                          << std::to_string(j) << ") = ";
                std::cout << M_gridEpetra(i,j)->norm1() << std::endl;
            }
        }
    }
}

} // namespace RedMA
