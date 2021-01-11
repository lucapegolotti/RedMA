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
#include <redma/problem/DataContainer.hpp>
#include <redma/problem/ComparisonFEMvsRB.hpp>

using namespace RedMA;

std::vector<std::pair<std::string, shp<DataContainer>>>
generateDatafiles(EPETRACOMM comm)
{
    std::vector<std::pair<std::string, shp<DataContainer>>> retVec;

    std::vector<double> podtol_field0;
    podtol_field0.push_back(2e-3);
    podtol_field0.push_back(1e-3);
    podtol_field0.push_back(5e-4);
    podtol_field0.push_back(1.5e-4);

    std::vector<double> podtol_field1;
    podtol_field1.push_back(1e-5);
    // podtol_field1.push_back(1e-4);

    std::vector<int> usePrimalSupremizers;
    usePrimalSupremizers.push_back(1);

    std::vector<int> useDualSupremizers;
    useDualSupremizers.push_back(1);

    std::vector<int> couplingnMax;
    couplingnMax.push_back(6);

    std::vector<int> useExtrapolation;
    useExtrapolation.push_back(0);

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
                            shp<DataContainer> data(new DataContainer());
                            data->setDatafile("datafiles/data");
                            data->setVerbose(false);
                            data->finalize();

                            data->setValueDouble("rb/online/basis/podtol_field0", pdt0);
                            data->setValueDouble("rb/online/basis/podtol_field1", pdt1);
                            data->setValueInt("rb/online/basis/useprimalsupremizers", prsup);
                            data->setValueInt("rb/online/basis/usedualsupremizers", dusup);
                            data->setValueInt("coupling/nMax", nMax);
                            data->setValueInt("time_discretization/use_extrapolation", extr);
                            data->setValueInt("rb/online/approximatenonlinearterm", 0);

                            // Chrono chrono;
                            // chrono.start();
                            // shp<ProblemRB> newProblemRB(new ProblemRB(data, comm));
                            // double setupTimeRB = chrono.diff();

                            std::string description;
                            description = "podtol_field0," + std::to_string(pdt0) + "\n";
                            description += "podtol_field1," + std::to_string(pdt1) + "\n";
                            description += "useprimalsupremizers," + std::to_string(prsup) + "\n";
                            description += "usedualsupremizers," + std::to_string(dusup) + "\n";
                            description += "couplingnMax," + std::to_string(nMax) + "\n";
                            description += "useextrapolation," + std::to_string(extr) + "\n";
                            description += "approximatenonlinearterm,0\n";
                            description += "numbernonlinearterms,0\n";
                            std::pair<std::string,shp<DataContainer>> newPair;
                            newPair.first = description;
                            newPair.second = data;
                            retVec.push_back(newPair);
                        }
                    }
                }
            }
        }
    }


    std::vector<int> nnterms;
    nnterms.push_back(10);
    nnterms.push_back(20);
    nnterms.push_back(40);
    nnterms.push_back(80);
    nnterms.push_back(120);

    for (auto nnt : nnterms)
    {
        shp<DataContainer> data(new DataContainer());
        data->setDatafile("datafiles/data");
        data->setVerbose(false);
        data->finalize();

        data->setValueDouble("rb/online/basis/podtol_field0", 1e-3);
        data->setValueDouble("rb/online/basis/podtol_field1", 1e-5);
        data->setValueInt("rb/online/basis/useprimalsupremizers", 1);
        data->setValueInt("rb/online/basis/usedualsupremizers", 1);
        data->setValueInt("coupling/nMax", 6);
        data->setValueInt("time_discretization/use_extrapolation", 0);
        data->setValueInt("rb/online/approximatenonlinearterm", 1);
        data->setValueInt("rb/online/numbernonlinearterms", nnt);

        // Chrono chrono;
        // chrono.start();
        // shp<ProblemRB> newProblemRB(new ProblemRB(data, comm));
        // double setupTimeRB = chrono.diff();

        std::string description;
        description = "podtol_field0," + std::to_string(1e-3) + "\n";
        description += "podtol_field1," + std::to_string(1e-5) + "\n";
        description += "useprimalsupremizers," + std::to_string(1) + "\n";
        description += "usedualsupremizers," + std::to_string(1) + "\n";
        description += "couplingnMax," + std::to_string(6) + "\n";
        description += "useextrapolation," + std::to_string(0) + "\n";
        description += "approximatenonlinearterm,1\n";
        description += "numbernonlinearterms," + std::to_string(nnt) + "\n";
        std::pair<std::string,shp<DataContainer>> newPair;
        newPair.first = description;
        newPair.second = data;
        retVec.push_back(newPair);
    }

    return retVec;
}

int main(int argc, char **argv)
{
    // #ifdef HAVE_MPI
    // MPI_Init (nullptr, nullptr);
    // EPETRACOMM comm (new Epetra_MpiComm(MPI_COMM_WORLD));
    // #else
    // EPETRACOMM comm(new Epetra_SerialComm());
    // #endif
    //
    // DataContainer data;
    // data.setDatafile("datafiles/data");
    // data.setVerbose(false);
    // std::string inletDirichlet = "weak";
    // data.setValueString("bc_conditions/inletdirichlet", inletDirichlet);
    // data.finalize();
    // Chrono chrono;
    // chrono.start();
    // shp<ProblemFEM> femProblem(new ProblemFEM(data, comm));
    // double setupTimeFEM = chrono.diff();
    //
    // ComparisonFEMvsRB comparison(data, comm);
    // comparison.setProblemFEM(femProblem);
    // double runTimeFEM;
    //
    // std::string outdir = "solution_fem_reference/";
    // if (fs::esists(outdir))
    // {
    //     std::ifstream infile;
    //
    //     comparison.loadFEMSolution(outdir);
    //
    //     infile.open(outdir + "femRunTime.txt");
    //
    //     std::string line;
    //     std::getline(infile,line);
    //
    //     runTimeFEM = std::atof(line.c_str());
    //     infile.close();
    // }
    // else
    // {
    //     comparison.runFEM();
    //     runTimeFEM = comparison.getTimeFem();
    //
    //     data.setValueString("exporter/outdir", outdir);
    //     comparison.exportFEM(1);
    //
    //     comparison.dumpFEMSolution(outdir);
    //
    //     std::ofstream outfile;
    //
    //     outfile.open(outdir + "femRunTime.txt", std::ios_base::app);
    //     outfile << runTimeFEM << "\n";
    //     outfile.close();
    // }
    //
    // std::vector<std::pair<std::string, SHP(DataContainer)>> datacs = generateDatafiles(comm);
    //
    // fs::create_directory("rb_solutions");
    //
    // unsigned int index = 0;
    // for (auto curdata : datacs)
    // {
    //     std::string curDir = "rb_solutions/sol" + std::to_string(index) + "/";
    //     fs::create_directory(curDir);
    //
    //     Chrono chrono;
    //     chrono.start();
    //     shp<ProblemRB> rbProblem(new ProblemRB(*curdata.second, comm));
    //     double setupTimeRB = chrono.diff();
    //
    //     comparison.setProblemRB(rbProblem);
    //     try
    //     {
    //         comparison.runRB();
    //         double runTimeRB = comparison.getTimeRB();
    //
    //         outdir = curDir + "solution/";
    //         rbProblem->getData().setValueString("exporter/outdir", outdir);
    //         data.setValueString("exporter/outdir", outdir);
    //         comparison.exportRB(4);
    //
    //         outdir = curDir + "error/";
    //         femProblem->getData().setValueString("exporter/outdir", outdir);
    //         comparison.exportError();
    //
    //         std::ofstream outfile;
    //
    //         outfile.open(curDir + "simulationData.txt", std::ios_base::app);
    //         outfile << curdata.first;
    //         outfile << "rbSetupTime" << setupTimeRB << "\n";
    //         outfile << "rbRunTime," << runTimeRB << "\n";
    //         outfile << "femSetupTime," << setupTimeFEM << "\n";
    //         outfile << "femRunTime," << runTimeFEM << "\n";
    //         outfile.close();
    //     } catch (Exception e)
    //     {
    //         std::cout << "Exception: probably solver has not converged" << std::endl << std::flush;
    //     }
    //     index++;
    // }

    return 0;
}
