#include "SteadySolver.hpp"

namespace RedMA {

SteadySolver::
SteadySolver(const DataContainer &data) :
    M_data(data),
    M_systemSolver(data){}

SteadySolver::
SteadySolver(const DataContainer &data, shp<FunProvider> funProvider) :
    M_data(data),
    M_systemSolver(data),
    M_funProvider(funProvider)
{
    initializeStatisticsFile();
}

void
SteadySolver::
initializeStatisticsFile()
{
    M_statisticsFile = M_data("time_discretization/solverstatistics","none");

    if (std::strcmp(M_statisticsFile.c_str(), "none"))
    {
        std::ofstream outfile;

        outfile.open(M_statisticsFile, std::ios_base::app);
        outfile << "solvetime,numiterations,precsetup\n";
        outfile.close();
    }
}

void
SteadySolver::
dumpSolverStatistics(std::vector<SolverStatistics> statistics) const
{
    if (std::strcmp(M_statisticsFile.c_str(),"none"))
    {
        std::ofstream outfile;
        outfile.open(M_statisticsFile, std::ios_base::app);
        for (auto entry : statistics)
        {
            outfile << std::to_string(entry.M_solveTime) << ","
            << std::to_string(entry.M_numIterations) << ","
            << std::to_string(entry.M_precSetupTime) << "\n";
        }
        outfile.close();
    }
}

void
SteadySolver::
setInitialGuess(const std::string &ICpath)
{
    M_initialGuess = this->M_funProvider->getZeroVector();

    if (fs::exists(ICpath) && (std::strcmp(ICpath.c_str(), "")))
    {
        std::map<unsigned int, std::vector<shp<aVector>>> IC_map = this->M_funProvider->importSolution(ICpath);
        unsigned int nPrimalBlocks = IC_map.size();

        for (unsigned int cnt = 0; cnt<nPrimalBlocks; cnt++)
            spcast<BlockVector>(M_initialGuess)->setBlock(cnt, IC_map[cnt][0]);
    }
}

void
SteadySolver::
setInitialGuess(const BV IC)
{
    M_initialGuess = IC;
}

shp<aVector>
SteadySolver::
solve(int& status)
{
    this->M_funProvider->applyDirichletBCs(0.0, M_initialGuess);

    FunctionFunctor<BV,BV> fct([this](BV sol)
    {
        BV retVec(this->M_funProvider->getRightHandSide(0.0, sol));

        // the previous solution satisfies the boundary conditions, so we search
        // for an increment with 0bcs
        this->M_funProvider->apply0DirichletBCs(retVec);

        return retVec;
    });

    FunctionFunctor<BV,BM> jac([this](BV sol)
    {
        // here the choice of hard copies is compulsory
        BM retMat(new BlockMatrix(1,1));
        retMat->deepCopy(this->M_funProvider->getJacobianRightHandSide(0.0, sol, 1.0));

        return retMat;
    });

    BV sol = this->M_systemSolver.solve(fct, jac, M_initialGuess, status);

    if (status != 0)
        throw new Exception("Solver has not converged!");

    std::vector<SolverStatistics> statistics = this->M_systemSolver.getSolverStatistics();
    this->dumpSolverStatistics(statistics);

    return sol;

}

}