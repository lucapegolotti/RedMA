namespace RedMA
{

template<class InVectorType, class InMatrixType>
aTimeMarchingAlgorithm<InVectorType,InMatrixType>::
aTimeMarchingAlgorithm(const DataContainer& data) :
  M_data(data),
  M_systemSolver(data)
{
}

template<class InVectorType, class InMatrixType>
aTimeMarchingAlgorithm<InVectorType,InMatrixType>::
aTimeMarchingAlgorithm(const DataContainer& data,
                       SHP(FunProvider) funProvider) :
  M_data(data),
  M_systemSolver(data),
  M_funProvider(funProvider)
{
    initializeStatisticsFile();
}

template<class InVectorType, class InMatrixType>
void
aTimeMarchingAlgorithm<InVectorType,InMatrixType>::
initializeStatisticsFile()
{
    M_statisticsFile = M_data("time_discretization/solverstatistics","none");

    if (std::strcmp(M_statisticsFile.c_str(),"none"))
    {
        std::ofstream outfile;

        outfile.open(M_statisticsFile, std::ios_base::app);
        outfile << "time,solvetime,numiterations,precsetup\n";
    }
}

template<class InVectorType, class InMatrixType>
void
aTimeMarchingAlgorithm<InVectorType,InMatrixType>::
dumpSolverStatistics(std::vector<SolverStatistics> statistics,
                     const double& t) const
{
    if (std::strcmp(M_statisticsFile.c_str(),"none"))
    {
        std::ofstream outfile;
        outfile.open(M_statisticsFile, std::ios_base::app);
        for (auto entry : statistics)
        {
            outfile << std::to_string(t) << ","
                    << std::to_string(entry.M_solveTime) << ","
                    << std::to_string(entry.M_numIterations) << ","
                    << std::to_string(entry.M_precSetupTime) << "\n";
        }
    }
}

}
