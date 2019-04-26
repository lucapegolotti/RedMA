// implementation of template class

namespace RedMA
{

template <class AssemblerType>
TimeMarchingAlgorithm<AssemblerType>::
TimeMarchingAlgorithm(const GetPot& datafile) :
  M_datafile(datafile)
{
}

}  // namespace RedMA
