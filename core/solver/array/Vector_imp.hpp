// implementation of Vector.hpp

namespace RedMA
{

template<class InVectorType>
Vector<InVectorType>::
Vector() :
  M_inVector(nullptr)
{
}

template<class InVectorType>
Vector<InVectorType>::
Vector(const InVectorTypePtr& matrix)
{
    M_inVector = matrix;
}

template<class InVectorType>
Vector<InVectorType>::
Vector(const Vector<InVectorType>& other)
{
    M_inVector = other.get();
}

template<class InVectorType>
Vector<InVectorType>&
Vector<InVectorType>::
operator=(const Vector<InVectorType>& other)
{
    M_inVector.reset(new InVectorType(*other.get()));
    return *this;
}

template<class InVectorType>
Vector<InVectorType>&
Vector<InVectorType>::
operator=(const InVectorTypePtr& vector)
{
    M_inVector.reset(new InVectorType(*vector));
    return *this;
}

template<class InVectorType>
Vector<InVectorType>&
Vector<InVectorType>::
operator+=(const Vector<InVectorType>& other)
{
    if (!isNull())
        *M_inVector += *other.get();
    else
        M_inVector.reset(new InVectorType(*other.get()));
    return *this;
}

template<class InVectorType>
Vector<InVectorType>&
Vector<InVectorType>::
operator*=(const double& coeff)
{
    if (!isNull())
        *M_inVector *= coeff;
    return *this;
}

template<class InVectorType>
typename Vector<InVectorType>::InVectorTypePtr&
Vector<InVectorType>::
get()
{
    return M_inVector;
}

template<class InVectorType>
typename Vector<InVectorType>::InVectorTypePtr
Vector<InVectorType>::
get() const
{
    return M_inVector;
}

template<class InVectorType>
bool
Vector<InVectorType>::
isZero() const
{
    return (M_inVector.norm1() < 5e-16);
}

template<class InVectorType>
bool
Vector<InVectorType>::
isNull() const
{
    return (M_inVector == nullptr);
}

}  // namespace RedMA
