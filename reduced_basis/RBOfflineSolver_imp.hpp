namespace RedMA
{

template<class AssemblerType>
RBOfflineSolver<AssemblerType>::
RBOfflineSolver(const GetPot& datafile, rbLifeV::ParameterHandler& parameterHandler,
                commPtr_Type comm, bool verbose) :
  GlobalSolver<AssemblerType>(datafile, comm, verbose),
  rbLifeV::ApproximatedAffineFemProblem(comm, parameterHandler)
{

}

}
