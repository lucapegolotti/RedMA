// implementation of Matrix.hpp

template <InMatrixPtr>
Matrix<InMatrixPtr>::
Matrix()
{
}

template <InMatrixPtr>
Matrix<InMatrixPtr>::
Matrix(unsigned int numRows, unsigned int numCols)
{
    resize(numRows, numCols);
}

template <InMatrixPtr>
unsigned int
Matrix<InMatrixPtr>::
getNumberRows()
{
    return M_rows;
}

template <InMatrixPtr>
unsigned int
Matrix<InMatrixPtr>::
getNumberCols()
{
    return M_cols;
}

template <InMatrixPtr>
void
Matrix<InMatrixPtr>::
zero()
{
    M_matrix = nullptr;
}

template <InMatrixPtr>
typename Matrix<InMatrixPtr>::InMatrixPtr
get()
{
    return M_matrix;
}
