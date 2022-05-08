//
// Created by Federico Betti on 4/6/22.
//

#include "SnapshotsSteadySampler.hpp"

namespace RedMA
{

SnapshotsSteadySampler::
SnapshotsSteadySampler(const DataContainer &data, EPETRACOMM comm, std::vector<unsigned int> numSamples) :
        M_data(data), M_comm(comm), M_StratifiedSampler(numSamples)
{
}

void
SnapshotsSteadySampler::
takeSnapshots(const unsigned int &Nstart)
{
    std::map <std::string, std::vector<double>> SnapshotsSamples =
            M_StratifiedSampler.generateSamples();

    std::string outdir = M_data("rb/offline/snapshots/directory", "snapshots");
    unsigned int numInletConditions = M_data("bc_conditions/numinletbcs", 1);
    unsigned int numOutletConditions = M_data("bc_conditions/numoutletbcs", 0);

    fs::create_directory(outdir);
    GeometryPrinter printer;

    double elapsedTime = 0.0;

    BV solutionPrevAmplitude;
    BV solutionPrevWidth;
    BV solutionPrevFlowRate;

    std::map<std::string, BV> initialConditions;

    for (unsigned int i = 0; i < M_StratifiedSampler.getNumComponents()[0]; i++)
    {
        std::map<std::string, double> currentSample;
        currentSample["flow_rate"] = SnapshotsSamples["flow_rate"][i];
        for (unsigned int j = 0; j < M_StratifiedSampler.getNumComponents()[1]; j++)
        {
            currentSample["stenosis_amplitude"] = SnapshotsSamples["stenosis_amplitude"][j];
            for (unsigned int k = 0; k < M_StratifiedSampler.getNumComponents()[2]; k++)
            {
                Chrono chrono;
                chrono.start();

                currentSample["stenosis_width"] = SnapshotsSamples["stenosis_width"][k];

                std::string msg = "Performing the snapshot generation for the current values";
                msg = msg + "The current parameters are: ";
                printlog(GREEN, msg);
                printCurrentSample(currentSample);

                std::vector<double> paramsVec = getParametersValuesAsVector(currentSample);

                // setting directory where to store the solutions
                std::string curdir = outdir + "/param" + std::to_string(i) + "_"
                                     + std::to_string(j) + "_" + std::to_string(k);

                // setting the flow rate at the unique inlet, the flow_rate is in position zero
                unsigned int numInlet = 0;
                std::function<double(double)> inletDirichlet = [paramsVec, numInlet](double t) {return paramsVec[numInlet]; };
                M_data.setInletBC(inletDirichlet, numInlet);

                GlobalProblem problem(M_data, M_comm, false);
                fs::create_directory(curdir);

                std::string curSolutionDir = curdir + "/";
                M_data.setValueString("exporter/outdir", curSolutionDir);

                problem.getTree().setGeometricParametersFromSample(currentSample);

                problem.setup();

                std::string filename = curdir + "/bypass.xml";
                printer.saveToFile(problem.getTree(), filename, M_comm);

                if (k != 0)
                    problem.getSteadySolver()->setInitialGuess(initialConditions["stenosis_width"]);
                if (j != 0 && k == 0)
                    problem.getSteadySolver()->setInitialGuess(initialConditions["stenosis_amplitude"]);
                if (i != 0 && j == 0 && k == 0)
                    problem.getSteadySolver()->setInitialGuess(initialConditions["flow_rate"]);
;
                problem.solveSteady();

                if (k == 0)
                    initialConditions["stenosis_amplitude"] = problem.getLastSolution();
                if (j == 0 && k == 0)
                    initialConditions["flow_rate"] = problem.getLastSolution();

                // this is always updated
                initialConditions["stenosis_width"] = problem.getLastSolution();

                elapsedTime += chrono.diff() / std::accumulate(begin(M_StratifiedSampler.getNumComponents()),
                                                               end(M_StratifiedSampler.getNumComponents()),
                                                               1, std::multiplies<unsigned int>());

                dumpSnapshots(problem, curSolutionDir, paramsVec);
            }
        }
    }

    std::string msg = "Average time per snapshot =  ";
    msg += std::to_string(elapsedTime);
    msg += " seconds\n";
    printlog(MAGENTA, msg, true);
}


//    for (unsigned int i = 0; i < M_StratifiedSampler.getNumSamples(); i++)
//    {
//        std::map<std::string, double> vec = SnapshotsSamples[i];
//        std::string msg = "Performing the # " + std::to_string(i) + "snapshot generation";
//        msg = msg + "The current parameters are: ";
//        printlog(GREEN, msg);
//        printCurrentSample(vec);
//
//        std::vector<double> paramsVec = getParametersValuesAsVector(vec);
//
//        // paramsVec.insert(paramsVec.begin() + 1, 1.0 - paramsVec[0]);
//
//        // to guarantee (almost...) that two snapshots are not saved at the same location!
//        unsigned int paramIndex = Nstart;
//        // we find the first parameter index available, starting from Nstart
//        while (fs::exists(outdir + "/param" + std::to_string(paramIndex)))
//            paramIndex++;
//        std::string curdir = outdir + "/param" + std::to_string(paramIndex);
//
//        for (unsigned int numInlet = 0; numInlet < numInletConditions; ++numInlet)
//        {
//            std::function<double(double)> inletDirichlet = [paramsVec, numInlet](double t) {return paramsVec[numInlet]; };
//            M_data.setInletBC(inletDirichlet, numInlet);
//        }
//
//       unsigned int numOutlet = 0;
//       M_data.setOutletBC([numOutlet](double t){return 1.5 * 1333.0;}, numOutlet);
//
//        // paramsVec.erase(paramsVec.begin() + 1);
//
//        GlobalProblem problem(M_data, M_comm, false);
//
//        problem.doStoreSolutions();
//
//        fs::create_directory(curdir);
//
//        std::string curSolutionDir = curdir + "/";
//        M_data.setValueString("exporter/outdir", curSolutionDir);
//
//        problem.getTree().setGeometricParametersFromSample(vec);
//
//        problem.setup();
//
//         if (!problem.isFEProblem())
//            throw new Exception("The tree must be composed of only FE nodes to "
//                                "sample the snapshots!");
//
//        std::string filename = curdir + "/bypass.xml";
//        printer.saveToFile(problem.getTree(), filename, M_comm);
//
//        problem.solveSteady();
//
//        dumpSnapshots(problem, curdir, getParametersValuesAsVector(vec));
//    }
//
//    std::string msg = "Average time per snapshot =  ";
//    msg += std::to_string(elapsedTime);
//    msg += " seconds\n";
//    printlog(MAGENTA, msg, true);
//
//}

void
SnapshotsSteadySampler::
dumpSnapshots(GlobalProblem& problem,
              std::string outdir, const std::vector<double>& params = {})
{
    auto IDmeshTypeMap = problem.getBlockAssembler()->getIDMeshTypeMap();
    auto solutions = problem.getSolutions();

    assert(!(solutions.empty()));

    unsigned int n_primal_blocks = IDmeshTypeMap.size();
    unsigned int n_dual_blocks = solutions[0]->nRows() - n_primal_blocks;

    auto M_mass = problem.getBlockAssembler()->block(0)->assembleMatrix(0);
    auto M_stiffness = problem.getBlockAssembler()->block(0)->assembleMatrix(1);
    auto M_divergence = problem.getBlockAssembler()->block(0)->assembleMatrix(2);

    unsigned int takeEvery = M_data("rb/offline/snapshots/take_every", 5);
    bool binary = M_data("rb/offline/snapshots/dumpbinary", true);
    bool computereynolds = M_data("rb/offline/snapshots/computereynolds", true);
    double density = M_data("fluid/density", 1.0);
    double viscosity = M_data("fluid/viscosity", 0.035);

    std::ios_base::openmode omode = std::ios_base::app;
    if (binary)
        omode = omode | std::ios::binary;

    M_mass->block(0, 0)->dump("M");
    M_stiffness->block(0, 0)->dump("A");
    M_divergence->block(0, 1)->dump("BdivT");
    M_divergence->block(1, 0)->dump("Bdiv");

    /*for (auto sol : solutions)
        problem.getBlockAssembler()->applyPiola(sol, true);*/

    for (auto idmeshtype: IDmeshTypeMap) {
        std::string meshtypedir = outdir + "/" + idmeshtype.second;
        fs::create_directory(meshtypedir);

        unsigned int nfields = solutions[0]->block(idmeshtype.first)->nRows();

        for (unsigned int i = 0; i < nfields; i++) {
            std::string outfilename = meshtypedir + "/field" + std::to_string(i) + ".snap";
            std::ofstream outfile;
            outfile.open(outfilename, omode);
            unsigned int count = 0;
            for (auto sol: solutions) {
                auto solBlck = convert<DistributedVector>(convert<BlockVector>(
                        convert<BlockVector>(sol)->block(idmeshtype.first))->block(i));
                if (count % takeEvery == 0) {
                    std::string str2write = solBlck->getString(',') + "\n";
                    if (M_comm->MyPID() == 0)
                        outfile.write(str2write.c_str(), str2write.size());
                }

                count++;
            }
            outfile.close();
        }

        for (unsigned int j = 0; j < n_dual_blocks; j++)  // save Lagrange multipliers, if any

        {
            std::string outfilename = meshtypedir + "/lagmult_" + std::to_string(j) + ".snap";
            std::ofstream outfile;
            outfile.open(outfilename, omode);
            unsigned int count = 0;
            unsigned int cur_index = n_primal_blocks + j;
            for (auto sol: solutions) {
                auto solBlck = convert<DistributedVector>(convert<BlockVector>(
                        convert<BlockVector>(sol)->block(cur_index))->block(0));
                if (count % takeEvery == 0) {
                    std::string str2write = solBlck->getString(',') + "\n";
                    if (M_comm->MyPID() == 0)
                        outfile.write(str2write.c_str(), str2write.size());
                }
                count++;
            }
            outfile.close();
        }

        shp<VECTOREPETRA> WSS(new VECTOREPETRA(problem.getBlockAssembler()->block(0)->getFESpaceBCs()->map, LifeV::Unique));
        if (computereynolds) {
            std::ofstream reynoldsfile;
            reynoldsfile.open(outdir + "/reynolds.txt", std::ios_base::app |
                                                        std::ios::binary);
            for (auto sol: solutions) {
                auto solBlck = convert<DistributedVector>(convert<BlockVector>(
                        convert<BlockVector>(sol)->block(idmeshtype.first))->block(0));
                double Umax = solBlck->maxMagnitude3D();
                auto tNode = problem.getBlockAssembler()->block(0)->getTreeNode();
                double D = 2 * tNode->M_block->getInlet(0).M_radius;
                double curReynolds = Umax * density * D / viscosity;

                if (M_comm->MyPID() == 0)
                    reynoldsfile << std::to_string(curReynolds) << std::endl;
            }
            reynoldsfile.close();
        }

        std::ofstream WSSfile;
        WSSfile.open(outdir + "/WSS.txt", std::ios_base::app |
                                          std::ios::binary);
        for (auto sol: solutions)
        {
            auto solBlck = convert<DistributedVector>(convert<BlockVector>(
                    convert<BlockVector>(sol)->block(idmeshtype.first))->block(0));
            problem.getBlockAssembler()->block(0)->computeWallShearStress(solBlck, WSS, M_comm);
            if (M_comm->MyPID() == 0)
                WSSfile << WSS << std::endl;
        }
        WSSfile.close();

        if (!params.empty()) {
            std::ofstream file(outdir + "/coeffile.txt", std::ios_base::app);

            for (auto element : params)
                file << std::fixed << std::setprecision(10) << element << std::endl;

            file.close();
        }
    }
}

void
SnapshotsSteadySampler::
printCurrentSample(std::map<std::string, double> sample)
{
    for (auto it = sample.begin(); it != sample.end(); ++it)
    {
        std::cout << it->first << " = " << it->second << std::endl;
    }
}

std::vector<double>
SnapshotsSteadySampler::
getParametersValuesAsVector(const std::map<std::string, double>& sample)
{

    std::vector<double> paramsValues;
    for (auto it = sample.begin(); it != sample.end(); ++ it)
    {
        paramsValues.push_back(it->second);
    }
    return paramsValues;
}
}