// implementation of Matrix.hpp

namespace RedMA
{

template<class InMatrixType>
Matrix<InMatrixType>::
Matrix() :
  M_inMatrix(nullptr)
{
}

template<class InMatrixType>
Matrix<InMatrixType>::
Matrix(const InMatrixTypePtr& matrix)
{
    M_inMatrix = matrix;
}

template<class InMatrixType>
Matrix<InMatrixType>::
Matrix(const Matrix<InMatrixType>& other)
{
    M_inMatrix = other.get();
}

template<class InMatrixType>
Matrix<InMatrixType>&
Matrix<InMatrixType>::
operator=(const Matrix<InMatrixType>& other)
{
    M_inMatrix.reset(new InMatrixType(*other.get()));
    return *this;
}

template<class InMatrixType>
Matrix<InMatrixType>&
Matrix<InMatrixType>::
operator=(const InMatrixTypePtr& matrix)
{
    M_inMatrix.reset(new InMatrixType(*matrix));
    return *this;
}

template<class InMatrixType>
Matrix<InMatrixType>&
Matrix<InMatrixType>::
operator+=(const Matrix<InMatrixType>& other)
{
    if (!isNull())
        *M_inMatrix += *other.get();
    else
        M_inMatrix.reset(new InMatrixType(*other.get()));
    return *this;
}

template<class InMatrixType>
Matrix<InMatrixType>&
Matrix<InMatrixType>::
operator*=(const double& coeff)
{
    if (!isNull())
        *M_inMatrix *= coeff;
    return *this;
}

template<class InMatrixType>
template<class InVectorType>
Vector<InVectorType>
Matrix<InMatrixType>::
operator*(const Vector<InVectorType>& vector)
{
    Vector<InVectorType> retVec;
    if (!isNull())
    {
        InVectorType res = (*M_inMatrix) * (*vector.get());
        retVec = std::make_shared<InVectorType>(new InVectorType(res));
    }
    return retVec;
}

template<class InMatrixType>
typename Matrix<InMatrixType>::InMatrixTypePtr&
Matrix<InMatrixType>::
get()
{
    return M_inMatrix;
}

template<class InMatrixType>
typename Matrix<InMatrixType>::InMatrixTypePtr
Matrix<InMatrixType>::
get() const
{
    return M_inMatrix;
}

template<class InMatrixType>
bool
Matrix<InMatrixType>::
isZero() const
{
    if (isNull())
        return false;
    return (M_inMatrix->norm1() < 5e-16);
}

template<class InMatrixType>
bool
Matrix<InMatrixType>::
isNull() const
{
    return (M_inMatrix == nullptr);
}

}  // namespace RedMA
