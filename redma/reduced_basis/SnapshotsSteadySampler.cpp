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
    std::map<std::string, std::vector<double>> SnapshotsSamples = M_StratifiedSampler.generateSamples();

    unsigned int branch_num = M_data("exporter/branch", 0);
    std::string outdir = M_data("rb/offline/snapshots/directory", "snapshots_" + std::to_string(branch_num));
    unsigned int numInletConditions = M_data("bc_conditions/numinletbcs", 1);
    unsigned int numOutletConditions = M_data("bc_conditions/numoutletbcs", 0);

    fs::create_directory(outdir);
    GeometryPrinter printer;

    double elapsedTime = 0.0;

    double DeltaT = 1.0;
    double curTime = 0.0;

    std::string solutionsDir = outdir + "/";
    M_data.setValueString("exporter/outdir", solutionsDir);
    M_data.setValueString("geometric_structure/xmlfile", "geometries/bypass" + std::to_string(branch_num) + ".xml");

    std::string assemblerType = M_data("assembler/type", "navierstokes");

    unsigned int numTotalSimulations = M_StratifiedSampler.getNumComponents()[0] *
                                       M_StratifiedSampler.getNumComponents()[1] *
                                       M_StratifiedSampler.getNumComponents()[2];

    // set correct parameters for the bypass
    tinyxml2::XMLDocument doc;
    int status = doc.LoadFile(("geometries/bypass" + std::to_string(branch_num) + ".xml").c_str());
    tinyxml2::XMLElement* element = doc.FirstChildElement();

    bool BL = false;
    if (element->Attribute("BL"))
    {
        BL = std::atoi(element->Attribute("BL"));
    }

    // returnBlock.reset(new Bypass(M_comm, "bypass", M_verbose, BL));

    bool isBifurcation = false;
    if (element->Attribute("isBifurcation"))
    {
        isBifurcation = std::atoi(element->Attribute("isBifurcation"));
    }

    // returnBlock.reset(new Bypass(M_comm, "bypass", M_verbose, BL, isBifurcation));

    unsigned int activeStenosis;
    if (element->Attribute("activeStenosis"))
    {
        activeStenosis = std::atoi(element->Attribute("activeStenosis"));
    }

    shp<Bypass> defaultBypass(new Bypass(M_comm, "bypass", false, BL, isBifurcation, activeStenosis, true));
    defaultBypass->readMesh();

    defaultBypass->setDiscretizationMethod("fem");
    defaultBypass->setAssemblerType(assemblerType);

    shp<TreeNode> defTreeNode(new TreeNode(defaultBypass, 0));
    shp<aAssembler> defAssembler = AssemblerFactory(M_data, defTreeNode);
    defAssembler->setup();

    std::map<std::string, double> currentSample;
    std::map<std::string, BV> initialConditions;
    BV TargetGuess;

    unsigned int count = 0;

    for (unsigned int i = 0; i < M_StratifiedSampler.getNumComponents()[0]; i++) {
        currentSample["flow_rate"] = SnapshotsSamples["flow_rate"][i];
        currentSample["stenosis_amplitude"] = 0;
        currentSample["stenosis_width"] = 0;

        std::vector<double> paramsVec = getParametersValuesAsVector(currentSample);
        unsigned int numInlet = 0;
        std::function<double(double)> inletDirichlet = [paramsVec, numInlet](
                double t) { return paramsVec[numInlet]; };
        M_data.setInletBC(inletDirichlet, numInlet);

        GlobalProblem problem(M_data, M_comm, false);
        problem.getTree().setGeometricParametersFromSample(currentSample);
        problem.setup();
        if (i != 0)
            problem.getSteadySolver()->setInitialGuess(TargetGuess);
        problem.solveSteady();
        TargetGuess = problem.getLastSolution();

        shp <BlockVector> solCopy(new BlockVector(2));
        // velocity
        solCopy->setBlock(0, problem.getLastSolution()->block(0)->block(0));
        // pressure
        solCopy->setBlock(1, problem.getLastSolution()->block(0)->block(1));

        std::map<unsigned int, double> targetRates = problem.getBlockAssembler()->block(0)->computeFlowRates(solCopy, true);

        for (auto it = targetRates.begin(); it != targetRates.end(); ++it)
        {
            std::cout << "Target rates : " << std::to_string(it->second) << endl;
        }

        for (unsigned int j = 0; j < M_StratifiedSampler.getNumComponents()[1]; j++) {
            currentSample["stenosis_amplitude"] = SnapshotsSamples["stenosis_amplitude"][j];
            for (unsigned int k = 0; k < M_StratifiedSampler.getNumComponents()[2]; k++) {

                Chrono chrono;
                chrono.start();

                currentSample["stenosis_width"] = SnapshotsSamples["stenosis_width"][k];

                std::string msg = "Performing the snapshot generation for the current values: \n";
                msg = msg + "The current parameters are: \n";
                printlog(GREEN, msg);
                printCurrentSample(currentSample);

                std::vector<double> paramsVec = getParametersValuesAsVector(currentSample);

                // setting directory where to store the solutions
                // std::string curdir = outdir + "/param" + std::to_string(i) + "_"
                //                     + std::to_string(j) + "_" + std::to_string(k);

                // setting the flow rate at the unique inlet, the flow_rate is in position zero
                unsigned int numInlet = 0;
                std::function<double(double)> inletDirichlet = [paramsVec, numInlet](
                        double t) { return paramsVec[numInlet]; };
                M_data.setInletBC(inletDirichlet, numInlet);

                // unsigned int numOutlet = 1;
                // std::function<double(double)> outletDirichlet = [paramsVec, numOutlet](
                //        double t) { return paramsVec[numOutlet]}; // moltiplicare per altro valore
                // M_data.setOutletBC(outletDirichlet, numOutlet);

                GlobalProblem problem(M_data, M_comm, false);
                // fs::create_directory(curdir);

                // std::string curSolutionDir = curdir + "/";

                problem.getTree().setGeometricParametersFromSample(currentSample);

                problem.doStoreSolutions();

                problem.setup();

                // std::string filename = curdir + "/bypass.xml";
                // printer.saveToFile(problem.getTree(), filename, M_comm);

                if (k != 0)
                    problem.getSteadySolver()->setInitialGuess(initialConditions["stenosis_width"]);
                if (j != 0 && k == 0)
                    problem.getSteadySolver()->setInitialGuess(initialConditions["stenosis_amplitude"]);
                if (i != 0 && j == 0 && k == 0)
                    problem.getSteadySolver()->setInitialGuess(initialConditions["flow_rate"]);

                problem.solveSteady();

                auto IDmeshTypeMap = problem.getBlockAssembler()->getIDMeshTypeMap();

                curTime += DeltaT;

                if (k == 0)
                    initialConditions["stenosis_amplitude"] = problem.getLastSolution();
                if (j == 0 && k == 0)
                    initialConditions["flow_rate"] = problem.getLastSolution();

                initialConditions["stenosis_width"] = problem.getLastSolution();

                shp <BlockVector> solCopy(new BlockVector(3));
                // velocity
                solCopy->setBlock(0, problem.getLastSolution()->block(0)->block(0));
                // pressure
                solCopy->setBlock(1, problem.getLastSolution()->block(0)->block(1));
                // displacement
                shp<VECTOREPETRA> disp = problem.getBlockAssembler()->block(0)->getTreeNode()->M_block->getDisplacement();
                shp<DistributedVector> dispVector(new DistributedVector());
                dispVector->setVector(disp);
                solCopy->setBlock(2, dispVector);

                defAssembler->exportSolution(curTime, solCopy);

                elapsedTime += chrono.diff() / numTotalSimulations;

                bool computereynolds = M_data("rb/offline/snapshots/computereynolds", true);
                double density = M_data("fluid/density", 1.0);
                double viscosity = M_data("fluid/viscosity", 0.035);

                for (auto idmeshtype: IDmeshTypeMap)
                {
                    if (computereynolds) {
                        std::ofstream reynoldsfile;
                        reynoldsfile.open(outdir + "/reynolds_" + std::to_string(branch_num) + ".txt", std::ios_base::app |
                                                                    std::ios::binary);
                        auto sol = problem.getLastSolution();
                        auto solBlck = convert<DistributedVector>(convert<BlockVector>(
                                convert<BlockVector>(sol)->block(idmeshtype.first))->block(0));
                        // problem.getBlockAssembler()->block(0)->getTreeNode()->computeWallShearStress(solBlck, WSS, M_comm);
                        double Umax = solBlck->maxMagnitude3D();
                        auto tNode = problem.getBlockAssembler()->block(0)->getTreeNode();
                        double D = 2 * tNode->M_block->getInlet(0).M_radius;
                        double curReynolds = Umax * density * D / viscosity;

                        if (M_comm->MyPID() == 0)
                            reynoldsfile << std::to_string(curReynolds) << std::endl;
                        reynoldsfile.close();
                    }
                }

                std::ofstream coeffile(outdir + "/coeffile_" + std::to_string(branch_num) + ".txt", std::ios_base::app);
                coeffile << std::fixed << std::setprecision(10) << "SIMULATION # " << std::to_string(count) << std::endl;
                coeffile << std::fixed << std::setprecision(10) << "Flow rate : " << std::to_string(currentSample["flow_rate"])
                     << std::endl;
                coeffile << std::fixed << std::setprecision(10) << "Stenosis amplitude : "
                     << std::to_string(currentSample["stenosis_amplitude"]) << std::endl;
                coeffile << std::fixed << std::setprecision(10) << "Stenosis width : "
                     << std::to_string(currentSample["stenosis_width"]) << std::endl;
                coeffile << std::endl << std::endl;
                coeffile.close();

                std::map<unsigned int, double> currRates = problem.getBlockAssembler()->block(0)->computeFlowRates(solCopy, false);

                std::ofstream ratesfile;
                ratesfile.open(outdir + "/flow_rates_" + std::to_string(branch_num) + ".txt", std::ios_base::app);
                ratesfile << std::fixed << std::setprecision(10) << "SIMULATION # "  << std::to_string(count) << std::endl;
                for (auto it = currRates.begin(); it != currRates.end(); ++it)
                {
                    ratesfile << std::fixed << std::setprecision(10) << std::to_string(it->first) << " : " << std::to_string(it->second) << std::endl;
                }
                ratesfile.close();

                std::ofstream paramsfile;
                paramsfile.open(outdir + "/parameters_" + std::to_string(branch_num) + ".txt", std::ios_base::app);
                paramsfile << std::fixed << std::setprecision(10) << std::to_string(currentSample["flow_rate"]) << ", ";
                paramsfile << std::fixed << std::setprecision(10) << std::to_string(currentSample["stenosis_amplitude"]) << ", ";
                paramsfile << std::fixed << std::setprecision(10) << std::to_string(currentSample["stenosis_width"]) << std::endl;
                paramsfile << std::endl;
                paramsfile.close();

                count += 1;

                std::ofstream FFRfile;
                FFRfile.open(outdir + "/FFR_" + std::to_string(branch_num) + ".txt", std::ios_base::app);
                for (auto it = currRates.begin(); it != currRates.end(); ++it)
                {
                    FFRfile << std::fixed << std::setprecision(10) << std::to_string(it->first) <<
                    " : " << std::to_string(it->second/targetRates[it->first]) << ", ";
                }
                FFRfile.close();
            }
        }
    }

    std::string msg = "Average time per snapshot =  ";
    msg += std::to_string(elapsedTime);
    msg += " seconds\n";
    printlog(MAGENTA, msg, true);
}

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

        if (computereynolds)
        {
            std::ofstream reynoldsfile;
            reynoldsfile.open(outdir + "/reynolds.txt", std::ios_base::app |
                                                        std::ios::binary);
            for (auto sol: solutions) {
                auto solBlck = convert<DistributedVector>(convert<BlockVector>(
                        convert<BlockVector>(sol)->block(idmeshtype.first))->block(0));
                // problem.getBlockAssembler()->block(0)->getTreeNode()->computeWallShearStress(solBlck, WSS, M_comm);
                double Umax = solBlck->maxMagnitude3D();
                auto tNode = problem.getBlockAssembler()->block(0)->getTreeNode();
                double D = 2 * tNode->M_block->getInlet(0).M_radius;
                double curReynolds = Umax * density * D / viscosity;

                if (M_comm->MyPID() == 0)
                    reynoldsfile << std::to_string(curReynolds) << std::endl;
            }
            reynoldsfile.close();
        }

        if (!params.empty()) {
            std::ofstream file(outdir + "/coeffile.txt", std::ios_base::app);

            for (auto element : params)
                file << std::fixed << std::setprecision(10) << element << std::endl;

            file.close();
        }
    }
}

//void
//SnapshotsSteadySampler::
//saveCoeffsFile(std::string outdir, std::map<std::string, std::vector<double>> samples)
//{
//    assert(!(samples.empty()));
//    unsigned int count = 0;
//
//    std::cout << "[SnapshotsSteadySampler] Saving all parameters to file..." << std::endl;
//    std::ofstream file(outdir + "/coeffile.txt", std::ios_base::app);
//
//    for (unsigned int i = 0; i < M_StratifiedSampler.getNumComponents()[0]; i++)
//    {
//        for (unsigned int j = 0; j < M_StratifiedSampler.getNumComponents()[1]; j++)
//        {
//            for (unsigned int k = 0; k < M_StratifiedSampler.getNumComponents()[2]; k++)
//            {
//                file << std::fixed << std::setprecision(10) << "SIMULATION # " << std::to_string(count) << std::endl;
//                file << std::fixed << std::setprecision(10) << "Flow rate : " << std::to_string(samples["flow_rate"][i]) << std::endl;
//                file << std::fixed << std::setprecision(10) << "Stenosis amplitude : " << std::to_string(samples["stenosis_amplitude"][j]) << std::endl;
//                file << std::fixed << std::setprecision(10) << "Stenosis width : " << std::to_string(samples["stenosis_width"][k]) << std::endl;
//                file << std::endl << std::endl;
//                count += 1;
//            }
//        }
//    }
//    file.close();
//}

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