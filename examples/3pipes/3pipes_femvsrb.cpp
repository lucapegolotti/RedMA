// Reduced Modeling of Arteries (RedMA)
// Copyright (C) 2019  Luca Pegolotti
//
// RedMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RedMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <redma/RedMA.hpp>
#include <redma/solver/problem/DataContainer.hpp>
#include <redma/solver/problem/ComparisonFEMvsRB.hpp>

using namespace RedMA;

std::vector<std::pair<std::string, SHP(ProblemRB)>>
generateRBproblems(DataContainer& data, EPETRACOMM comm)
{
    std::vector<std::pair<std::string, SHP(ProblemRB)>> retVec;

    std::vector<double> podtol_field0;
    podtol_field0.push_back(2e-3);
    podtol_field0.push_back(1e-3);
    podtol_field0.push_back(5e-4);
    podtol_field0.push_back(1e-4);

    std::vector<double> podtol_field1;
    podtol_field1.push_back(1e-5);
    // podtol_field1.push_back(1e-4);

    std::vector<std::string> usePrimalSupremizers;
    usePrimalSupremizers.push_back("true");

    std::vector<std::string> useDualSupremizers;
    useDualSupremizers.push_back("true");

    std::vector<int> couplingnMax;
    couplingnMax.push_back(6);

    std::vector<std::string> useExtrapolation;
    useExtrapolation.push_back("false");

    for (auto nMax : couplingnMax)
    {
        for (auto extr : useExtrapolation)
        {
            for (auto dusup : useDualSupremizers)
            {
                for (auto prsup : usePrimalSupremizers)
                {
                    for (auto pdt1 : podtol_field1)
                    {
                        for (auto pdt0 : podtol_field0)
                        {
                            data.setValue("rb/online/basis/podtol_field0", pdt0);
                            data.setValue("rb/online/basis/podtol_field1", pdt1);
                            data.setValue("rb/online/basis/useprimalsupremizers", prsup);
                            data.setValue("rb/online/basis/usedualsupremizers", dusup);
                            data.setValue("coupling/nMax", nMax);
                            data.setValue("time_discretization/use_extrapolation", extr);

                            LifeV::LifeChrono chrono;
                            chrono.start();
                            SHP(ProblemRB) newProblemRB(new ProblemRB(data, comm));
                            double setupTimeRB = chrono.diff();

                            std::string description;
                            description = "podtol_field0," + std::to_string(pdt0) + "\n";
                            description += "podtol_field1," + std::to_string(pdt1) + "\n";
                            description += "useprimalsupremizers," + std::to_string(prsup == "true") + "\n";
                            description += "usedualsupremizers," + std::to_string(dusup == "true") + "\n";
                            description += "couplingnMax," + std::to_string(nMax) + "\n";
                            description += "useextrapolation," + std::to_string(extr == "true") + "\n";
                            description += "rbSetupTime," + std::to_string(setupTimeRB) + "\n";
                            std::pair<std::string,SHP(ProblemRB)> newPair;
                            newPair.first = description;
                            newPair.second = newProblemRB;
                            retVec.push_back(newPair);
                        }
                    }
                }
            }
        }
    }

    return retVec;
}

int main(int argc, char **argv)
{
    #ifdef HAVE_MPI
    MPI_Init (nullptr, nullptr);
    EPETRACOMM comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    #else
    EPETRACOMM comm(new Epetra_SerialComm());
    #endif

    DataContainer data;
    data.setDatafile("datafiles/data");
    data.setVerbose(comm->MyPID() == 0);
    std::string inletDirichlet = "strong";
    data.setValue("bc_conditions/inletdirichlet", inletDirichlet);
    data.finalize();
    LifeV::LifeChrono chrono;
    chrono.start();
    SHP(ProblemFEM) femProblem(new ProblemFEM(data, comm));
    double setupTimeFEM = chrono.diff();

    ComparisonFEMvsRB comparison(data, comm);
    comparison.setProblemFEM(femProblem);
    double runTimeFEM;

    std::string outdir = "solution_fem_reference/";
    if (boost::filesystem::exists(outdir))
    {
        std::ifstream infile;

        comparison.loadFEMSolution(outdir);

        infile.open(outdir + "femRunTime.txt");

        std::string line;
        std::getline(infile,line);

        runTimeFEM = std::atof(line.c_str());
        infile.close();
    }
    else
    {
        comparison.runFEM();
        runTimeFEM = comparison.getTimeFem();

        data.setValue("exporter/outdir", outdir);
        comparison.exportFEM(1);

        comparison.dumpFEMSolution(outdir);

        std::ofstream outfile;

        outfile.open(outdir + "femRunTime.txt", std::ios_base::app);
        outfile << runTimeFEM << "\n";
        outfile.close();
    }

    // reset weak bcs for reduced basis
    inletDirichlet = "weak";
    data.setValue("bc_conditions/inletdirichlet", inletDirichlet);

    std::vector<std::pair<std::string, SHP(ProblemRB)>> rbProblems = generateRBproblems(data, comm);

    boost::filesystem::create_directory("rb_solutions");

    unsigned int index = 0;
    for (auto rbProblem : rbProblems)
    {
        std::string curDir = "rb_solutions/sol" + std::to_string(index) + "/";
        boost::filesystem::create_directory(curDir);

        comparison.setProblemRB(rbProblem.second);
        try
        {
            comparison.runRB();
            double runTimeRB = comparison.getTimeRB();

            outdir = curDir + "solution/";
            data.setValue("exporter/outdir", outdir);
            comparison.exportRB(4);

            outdir = curDir + "error/";
            data.setValue("exporter/outdir", outdir);
            comparison.exportError();

            std::ofstream outfile;

            outfile.open(curDir + "simulationData.txt", std::ios_base::app);
            outfile << rbProblem.first;
            outfile << "rbRunTime," << runTimeRB << "\n";
            outfile << "femSetupTime," << setupTimeFEM << "\n";
            outfile << "femRunTime," << runTimeFEM << "\n";
            outfile.close();
        } catch (Exception e)
        {
            std::cout << "Exception: probably solver has not converged" << std::endl << std::flush;
        }
        index++;
    }

    return 0;
}
