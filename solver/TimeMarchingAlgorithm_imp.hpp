// implementation of template class

namespace RedMA
{

template <class AssemblerType>
TimeMarchingAlgorithm<AssemblerType>::
TimeMarchingAlgorithm(const GetPot& datafile,
                      GlobalAssemblerType* assembler) :
  M_datafile(datafile),
  M_globalAssembler(assembler)
{
    M_solution.reset(new Vector(assembler->getGlobalMap()));
    M_solution->zero();
}

}  // namespace RedMA
