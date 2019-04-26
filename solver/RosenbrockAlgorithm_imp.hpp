// implementation of template class

namespace RedMA
{

template <class AssemblerType>
RosenbrockAlgorithm<AssemblerType>::
RosenbrockAlgorithm(const GetPot& datafile) :
  TimeMarchingAlgorithm<AssemblerType>(datafile)
{
}

template <class AssemblerType>
void
RosenbrockAlgorithm<AssemblerType>::
solveTimestep(const double &time, double &dt,
              const GlobalAssemblerType& assembler)
{
    std::cout << "Hey I am solving!" << std::endl;
}

}  // namespace RedMA
