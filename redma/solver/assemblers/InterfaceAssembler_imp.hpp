namespace RedMA
{


template <class InVectorType, class InMatrixType>
Interface<InVectorType, InMatrixType>::
Interface(SHP(AssemblerType) assemblerFather, const unsigned int& indexFather,
          SHP(AssemblerType) assemblerChild, const unsigned int& indexChild) :
  M_assemblerFather(assemblerFather),
  M_indexFather(indexFather),
  M_assemblerChild(assemblerChild),
  M_indexChild(indexChild)
{

}
}
