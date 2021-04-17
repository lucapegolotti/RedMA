#include "ComparisonFEMvsRB.hpp"

namespace RedMA
{

ComparisonFEMvsRB::
ComparisonFEMvsRB(const DataContainer& data, EPETRACOMM comm) :
  M_comm(comm),
  M_data(data)
{

}

void
ComparisonFEMvsRB::
runFEM()
{
    // Chrono chrono;
    // chrono.start();
    //
    // std::string msg = "Starting chrono for fem problem\n";
    // printlog(MAGENTA, msg, true);
    //
    // // deactivate exporter
    // M_data.setValueInt("exporter/save_every", -1);
    //
    // M_problemFEM->doStoreSolutions();
    // M_problemFEM->solve();
    //
    // M_timeFEM = chrono.diff();
    // msg = "Total time =  ";
    // msg += std::to_string(M_timeFEM);
    // msg += " seconds\n";
    // printlog(MAGENTA, msg, true);
}

void
ComparisonFEMvsRB::
runRB()
{
    // Chrono chrono;
    // chrono.start();
    //
    // std::string msg = "Starting chrono for rb problem\n";
    // printlog(MAGENTA, msg, true);
    //
    // // deactivate exporter
    // M_data.setValueInt("exporter/save_every", -1);
    //
    // M_problemRB->doStoreSolutions();
    // M_problemRB->solve();
    //
    // M_timeRB = chrono.diff();
    // msg = "Total time =  ";
    // msg += std::to_string(M_timeRB);
    // msg += " seconds\n";
    // printlog(MAGENTA, msg, true);
}



void
ComparisonFEMvsRB::
exportError()
{
    // std::string msg = "[ComparisonFEMvsRB] Exporting error solution\n";
    // printlog(MAGENTA, msg, true);
    //
    // unsigned int numPrimalBlocks = M_problemRB->getBlockAssembler()->getAssemblersMap().size();
    //
    // // we do this to reset the output directory
    // M_problemFEM->getBlockAssembler()->setExporter();
    //
    // unsigned int nsolutions;
    // if (M_loadedSolutions.size() == 0)
    //     nsolutions = M_problemFEM->getSolutions().size();
    // else
    //     nsolutions = M_loadedSolutions.size();
    // for (unsigned int i = 0; i < nsolutions; i++)
    // {
    //     BlockVector<BlockVector<VectorEp>> diff;
    //     diff.shallowCopy(M_problemRB->getBlockAssembler()->convertFunctionRBtoFEM(M_problemRB->getSolutions()[i], M_comm));
    //
    //     // we do this so that we don't have problems with different number of lagrange multipliers
    //     for (unsigned int j = 0; j < numPrimalBlocks; j++)
    //     {
    //         if (M_loadedSolutions.size() == 0)
    //             diff.block(j) -= M_problemFEM->getSolutions()[i].block(j);
    //         else
    //             diff.block(j) -= M_loadedSolutions[i].block(j);
    //     }
    //
    //     double t = M_problemRB->getTimesteps()[i];
    //     M_problemFEM->getBlockAssembler()->exportSolution(t, diff);
    // }
}

void
ComparisonFEMvsRB::
exportFEM(unsigned int saveEvery)
{
    // std::string msg = "[ComparisonFEMvsRB] Exporting FEM solution\n";
    // printlog(MAGENTA, msg, true);
    //
    // double t0 = M_data("time_discretization/t0", 0.0);
    // double T = M_data("time_discretization/T", 1.0);
    // double dt = M_data("time_discretization/dt", 0.01);
    //
    // double t = t0;
    //
    // // we do this to reset the output directory
    // M_problemFEM->getBlockAssembler()->setExporter();
    // unsigned int nsolutions = M_problemFEM->getSolutions().size();
    // for (unsigned int i = 0; i < nsolutions; i++)
    // {
    //     if (i % saveEvery == 0)
    //     {
    //         auto solution = M_problemFEM->getSolutions()[i];
    //         double t = M_problemFEM->getTimesteps()[i];
    //         M_problemFEM->getBlockAssembler()->exportSolution(t, solution);
    //     }
    // }
}

void
ComparisonFEMvsRB::
dumpFEMSolution(std::string outdir)
{
    // fs::create_directory(outdir);
    //
    // auto solutions = M_problemFEM->getSolutions();
    // unsigned int nsolutions = solutions.size();
    // unsigned int numPrimalBlocks = M_problemFEM->getBlockAssembler()->getAssemblersMap().size();
    //
    // for (unsigned int iblock = 0; iblock < numPrimalBlocks; iblock++)
    // {
    //     std::string curDir = outdir + "/block" + std::to_string(iblock) + "/";
    //     fs::create_directory(outdir + "/block" + std::to_string(iblock));
    //     for (unsigned int ifield = 0; ifield < solutions[0].block(iblock).nRows(); ifield++)
    //     {
    //         std::ofstream outfile;
    //         outfile.open(curDir + "field" + std::to_string(ifield) + ".txt", std::ios_base::app);
    //
    //         for (unsigned int isol = 0; isol < nsolutions; isol++)
    //         {
    //             std::string str2write = solutions[isol].block(iblock).block(ifield).getString(',') + "\n";
    //             outfile.write(str2write.c_str(), str2write.size());
    //         }
    //
    //         outfile.close();
    //     }
    // }
}


void
ComparisonFEMvsRB::
loadFEMSolution(std::string indir)
{
    // auto assMap = M_problemFEM->getBlockAssembler()->getAssemblersMap();
    // unsigned int numPrimalBlocks = assMap.size();
    // unsigned int numFields = assMap[0]->getNumComponents();
    //
    // int numSolutions = -1;
    //
    // std::vector<std::vector<std::vector<shp<VECTOREPETRA>>>> allSols;
    //
    // for (unsigned int iblock = 0; iblock < numPrimalBlocks; iblock++)
    // {
    //     std::vector<std::vector<shp<VECTOREPETRA>>> stackSol;
    //     for (unsigned int ifield = 0; ifield < numFields; ifield++)
    //     {
    //         std::vector<shp<VECTOREPETRA>> curSols;
    //         std::ifstream infile(indir + "/block" + std::to_string(iblock) + "/field" + std::to_string(ifield) + ".txt");
    //         std::string line;
    //         unsigned int count = 0;
    //         while(std::getline(infile,line))
    //         {
    //             shp<VECTOREPETRA> newVector(new VECTOREPETRA(assMap[iblock]->getFEspace(ifield)->map()));
    //
    //             std::stringstream linestream(line);
    //             std::string value;
    //             unsigned int i = 0;
    //             while(getline(linestream,value,','))
    //             {
    //                newVector->operator[](i) = std::atof(value.c_str());
    //                i++;
    //             }
    //             if (i != newVector->epetraVector().GlobalLength())
    //                throw new Exception("Stored solution length does not match fespace dimension!");
    //
    //             curSols.push_back(newVector);
    //             count++;
    //         }
    //         if (numSolutions == -1)
    //             numSolutions = count;
    //         else if (count != numSolutions)
    //             throw new Exception("Inconsistent number of solutions!");
    //         infile.close();
    //         stackSol.push_back(curSols);
    //     }
    //     allSols.push_back(stackSol);
    // }
    //
    // M_loadedSolutions.resize(numSolutions);
    // for (unsigned int isol = 0; isol < numSolutions; isol++)
    // {
    //     M_loadedSolutions[isol].resize(numPrimalBlocks);
    //     for (unsigned int iblock = 0; iblock < numPrimalBlocks; iblock++)
    //     {
    //         M_loadedSolutions[isol].block(iblock).resize(numFields);
    //         for (unsigned int ifield = 0; ifield < numFields; ifield++)
    //         {
    //             M_loadedSolutions[isol].block(iblock).block(ifield).data() = allSols[iblock][ifield][isol];
    //         }
    //     }
    // }
}

void
ComparisonFEMvsRB::
exportRB(unsigned int saveEvery)
{
    // std::string msg = "[ComparisonFEMvsRB] Exporting RB solution\n";
    // printlog(MAGENTA, msg, true);
    //
    // // we do this to reset the output directory
    // M_problemRB->getBlockAssembler()->setExporter();
    //
    // unsigned int nsolutions = M_problemRB->getSolutions().size();
    // for (unsigned int i = 0; i < nsolutions; i++)
    // {
    //     if (i % saveEvery == 0)
    //     {
    //         auto solution = M_problemRB->getSolutions()[i];
    //         double t = M_problemRB->getTimesteps()[i];
    //         M_problemRB->getBlockAssembler()->exportSolution(t, solution);
    //     }
    // }
}

}
