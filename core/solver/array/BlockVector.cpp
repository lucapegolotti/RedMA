#include <BlockVector.hpp>

namespace RedMA
{

BlockVector::
BlockVector()
{
}

BlockVector::
BlockVector(unsigned int numRows)
{
    resize(numRows);
}

void
BlockVector::
resize(unsigned int numRows)
{
    M_vectors.resize(numRows);

    for (unsigned int i = 0; i < numRows; i++)
        M_vectors[i] = nullptr;

}

Vector<BlockVector::VectorEpetra>&
BlockVector::
block(unsigned int index)
{
    return M_vectors[index];
}

Vector<BlockVector::VectorEpetra>
BlockVector::
block(unsigned int index) const
{
    return M_vectors[index];
}

BlockVector::VectorEpetraPtr&
BlockVector::
operator[](unsigned int index)
{
    return block(index).get();
}

void
BlockVector::
push_back(Vector<VectorEpetra> newElement)
{
    M_vectors.push_back(newElement);
}

unsigned int
BlockVector::
size() const
{
    return M_vectors.size();
}


}  // namespace RedMA
