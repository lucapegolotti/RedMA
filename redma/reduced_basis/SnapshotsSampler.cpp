#include "SnapshotsSampler.hpp"

namespace RedMA

{

SnapshotsSampler::
SnapshotsSampler(const DataContainer& data, EPETRACOMM comm) :
  M_data(data),
  M_comm(comm)
{
}

void
SnapshotsSampler::
takeSnapshots(const unsigned int& Nstart)
{
    std::string outdir = M_data("rb/offline/snapshots/directory", "snapshots");
    std::string param_type = M_data("rb/offline/snapshots/param_type", "geometric");
    std::list<std::string> param_types = M_data.stringTokenizer(param_type, ',');

    fs::create_directory(outdir);
    GeometryPrinter printer;

    unsigned int nSnapshots = M_data("rb/offline/snapshots/number", 10);
    double bound = M_data("rb/offline/snapshots/bound", 0.2);

    double elapsedTime = 0.0;
    double elapsedTimeNoSetup = 0.0;

    for (unsigned int i = 0; i < nSnapshots; i++)
    {
        Chrono chrono;
        chrono.start();

        std::vector<double> array_params;

        // to guarantee (almost...) that two snapshots are not saved at the same location!
        unsigned int paramIndex = Nstart;
        
        // we find the first parameter index available, starting from Nstart
        while (fs::exists(outdir + "/param" + std::to_string(paramIndex)))
            paramIndex++;
        std::string curdir = outdir + "/param" + std::to_string(paramIndex);

        if (std::find(std::begin(param_types), std::end(param_types), "inflow") != std::end(param_types))
        {
            std::vector<double> array_params_inflow = this->sampleParametersInflow();
            array_params.insert(array_params.end(), array_params_inflow.begin(), array_params_inflow.end());
        }


        if (std::find(std::begin(param_types), std::end(param_types), "physics") != std::end(param_types))
        {
            std::vector<double> array_params_physics = this->sampleParametersPhysics();
            array_params.insert(array_params.end(), array_params_physics.begin(), array_params_physics.end());
        }

        GlobalProblem problem(M_data, M_comm, false);

        problem.doStoreSolutions();

        fs::create_directory(curdir);

        if (std::find(std::begin(param_types), std::end(param_types), "geometric") != std::end(param_types))
            problem.getTree().randomSampleAroundOriginalValue(bound);

        Chrono chrono2;
        chrono2.start();
        problem.setup();
        double setupTime = chrono2.diff();

        if (!problem.isFEProblem())
            throw new Exception("The tree must be composed of only FE nodes to "
                                "sample the snapshots!");

        std::string filename = curdir + "/tree.xml";
        printer.saveToFile(problem.getTree(), filename, M_comm);

        problem.solve();

        elapsedTime += chrono.diff() / nSnapshots;
        elapsedTimeNoSetup += (chrono.diff() - setupTime) / nSnapshots;

        dumpSnapshots(problem, curdir, array_params);
    }

    std::string msg = "Average time per snapshot =  ";
    msg += std::to_string(elapsedTime);
    msg += " seconds\n";
    printlog(MAGENTA, msg, true);

    msg = "Average time per snapshot (no setup) =  ";
    msg += std::to_string(elapsedTimeNoSetup);
    msg += " seconds\n";
    printlog(MAGENTA, msg, true);

}

void
SnapshotsSampler::
dumpSnapshots(GlobalProblem& problem,
              std::string outdir,
              const std::vector<double> array_params = {})
{
    auto IDmeshTypeMap = problem.getBlockAssembler()->getIDMeshTypeMap();

    auto solutions = problem.getSolutions();
    auto initialConditions = problem.getInitialConditions();
    auto extraSolutions = problem.getExtraSolutions();
    auto extraInitialConditions = problem.getExtraInitialConditions();

    assert(!(solutions.empty()));

    unsigned int n_primal_blocks = IDmeshTypeMap.size();
    unsigned int n_dual_blocks = solutions[0]->nRows() - n_primal_blocks;

    /*auto M_mass = problem.getBlockAssembler()->block(0)->assembleMatrix(0);
    auto M_stiffness = problem.getBlockAssembler()->block(0)->assembleMatrix(1);
    auto M_divergence = problem.getBlockAssembler()->block(0)->assembleMatrix(2);*/

    std::string param_type = M_data("rb/offline/snapshots/param_type", "geometric");
    std::list<std::string> param_types = M_data.stringTokenizer(param_type, ',');
    unsigned int takeEvery = M_data("rb/offline/snapshots/take_every", 5);
    bool binary = M_data("rb/offline/snapshots/dumpbinary", true);
    bool computereynolds = M_data("rb/offline/snapshots/computereynolds", true);
    double density = M_data("rb/offline/fluid/density", 1.06);
    double viscosity = M_data("rb/offline/fluid/viscosity", 0.035);

    std::ios_base::openmode omode = std::ios_base::app;
    if (binary)
        omode = omode | std::ios::binary;

    /*M_mass->block(0,0)->dump("M");
    M_stiffness->block(0,0)->dump("A");
    M_divergence->block(0,1)->dump("BdivT");
    M_divergence->block(1,0)->dump("Bdiv");*/

    if (std::find(std::begin(param_types), std::end(param_types), "geometric") != std::end(param_types))
        for (auto sol : solutions)
            problem.getBlockAssembler()->applyPiola(sol, true);

    for (auto idmeshtype : IDmeshTypeMap)
    {
        std::string meshtypedir = outdir + "/" + idmeshtype.second;
        fs::create_directory(meshtypedir);

        unsigned int nfields = solutions[0]->block(idmeshtype.first)->nRows();
        unsigned int nfields_extra = 0;
        if (extraSolutions.size())
            nfields_extra = extraSolutions[0]->block(idmeshtype.first)->nRows();

        for (unsigned int i = 0; i < nfields; i++)
        {
            std::string outfilename = meshtypedir + "/field" + std::to_string(i) + ".snap";

            std::ofstream outfile;
            outfile.open(outfilename, omode);
            unsigned int count = 0;
            for (auto sol : solutions)
            {
                auto solBlck = convert<DistributedVector>(convert<BlockVector>(
                               convert<BlockVector>(sol)->block(idmeshtype.first))->block(i));
                if (count % takeEvery == 0)
                {
                    std::string str2write = solBlck->getString(',') + "\n";
                    if (M_comm->MyPID() == 0)
                        outfile.write(str2write.c_str(), str2write.size());
                }

                count++;
            }
            outfile.close();
        }

        if (!(std::strcmp(M_data("assembler/type", "navierstokes").c_str(), "navierstokes_membrane")))
        {
            for(unsigned int k = 0; k < nfields_extra; k++)
                {
                std::string outfilename = meshtypedir + "/field" + std::to_string(nfields+k) + ".snap";
                std::ofstream outfile;
                outfile.open(outfilename, omode);
                unsigned int count = 0;
                for (auto sol : extraSolutions)
                {
                    auto solBlck = convert<DistributedVector>(convert<BlockVector>(
                            convert<BlockVector>(sol)->block(idmeshtype.first))->block(k));
                    if (count % takeEvery == 0)
                    {
                        std::string str2write = solBlck->getString(',') + "\n";
                        if (M_comm->MyPID() == 0)
                            outfile.write(str2write.c_str(), str2write.size());
                    }
                    count++;
                }
                outfile.close();
                }
        }

        if (!initialConditions.empty())
        {
            for (unsigned int i = 0; i < nfields; i++)
            {
                std::string outfilename = meshtypedir + "/field" + std::to_string(i) + "_IC.snap";

                std::ofstream outfile;
                outfile.open(outfilename, omode);
                unsigned int count = 0;
                for (auto initialCondition : initialConditions)
                {
                    auto solBlck = convert<DistributedVector>(convert<BlockVector>(
                            convert<BlockVector>(initialCondition)->block(idmeshtype.first))->block(i));
                    if (count % takeEvery == 0)
                    {
                        std::string str2write = solBlck->getString(',') + "\n";
                        if (M_comm->MyPID() == 0)
                            outfile.write(str2write.c_str(), str2write.size());
                    }

                    count++;
                }
                outfile.close();
            }
        }

        if (!extraInitialConditions.empty())
        {
            for (unsigned int k = 0; k < nfields_extra; k++)
            {
                std::string outfilename = meshtypedir + "/field" + std::to_string(nfields+k) + "_IC.snap";

                std::ofstream outfile;
                outfile.open(outfilename, omode);
                unsigned int count = 0;
                for (auto initialCondition : extraInitialConditions)
                {
                    auto solBlck = convert<DistributedVector>(convert<BlockVector>(
                            convert<BlockVector>(initialCondition)->block(idmeshtype.first))->block(k));
                    if (count % takeEvery == 0)
                    {
                        std::string str2write = solBlck->getString(',') + "\n";
                        if (M_comm->MyPID() == 0)
                            outfile.write(str2write.c_str(), str2write.size());
                    }

                    count++;
                }
                outfile.close();
            }
        }

        if (computereynolds)
        {
            std::ofstream reynoldsfile;
            reynoldsfile.open(outdir + "/reynolds.txt", std::ios_base::app |
                                                             std::ios::binary);
            for (auto sol : solutions)
            {
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

        if (!array_params.empty())
        {
            std::ofstream file(outdir + "/coeffile.txt", std::ios_base::app);

            for (double i : array_params)
                file << std::fixed << std::setprecision(10) << i << std::endl;

            file.close();
        }
    }

    // save Lagrange multipliers, if any
    std::string lagmultdir = outdir + "/multipliers";
    fs::create_directory(lagmultdir);
    for(unsigned int j = 0; j < n_dual_blocks; j++)
        {
        std::string outfilename = lagmultdir + "/lagmult" + std::to_string(j) + ".snap";
        std::ofstream outfile;
        outfile.open(outfilename, omode);
        unsigned int count = 0;
        unsigned int cur_index = n_primal_blocks + j;
        for (auto sol : solutions)
        {
            auto solBlck = convert<DistributedVector>(convert<BlockVector>(
                    convert<BlockVector>(sol)->block(cur_index))->block(0));
            if (count % takeEvery == 0)
            {
                std::string str2write = solBlck->getString(',') + "\n";
                if (M_comm->MyPID() == 0)
                    outfile.write(str2write.c_str(), str2write.size());
            }
            count++;
        }
        outfile.close();
        }
}

void
SnapshotsSampler::
transformSnapshotsWithPiola(std::string snapshotsDir,
                            unsigned int fieldIndex,
                            unsigned int maxSnapshot)
{
    std::string param_type = M_data("rb/offline/snapshots/param_type", "geometric");
    std::list<std::string> param_types = M_data.stringTokenizer(param_type, ',');
    if (std::find(std::begin(param_types), std::end(param_types), "geometric") == std::end(param_types))
        return;

    std::ios_base::openmode omode = std::ios_base::app;

    for (unsigned int i = 0; i < maxSnapshot; i++)
    {
        std::string paramfile = snapshotsDir + "/param" + std::to_string(i);
        std::string treefile =  paramfile + "/tree.xml";
        if (fs::exists(treefile))
        {
            M_data.setValueString("geometric_structure/xmlfile", treefile);
            GlobalProblem problem(M_data, M_comm, false);

            problem.setup();

            auto IDmeshTypeMap = problem.getBlockAssembler()->getIDMeshTypeMap();

            std::map<std::string, std::vector<unsigned int>> meshTypeNumber;

            for (auto m : IDmeshTypeMap)
                meshTypeNumber[m.second].push_back(m.first);

            for (auto mtn : meshTypeNumber)
            {
                auto auxVec = problem.getBlockAssembler()->getZeroVector();

                std::string snapshotFile = paramfile + "/" + mtn.first +
                                           "/field" + std::to_string(fieldIndex) + ".snap";
                std::cout << "snapshotFile = " << snapshotFile << std::endl << std::flush;
                std::ifstream infile(snapshotFile);
                std::string line;

                auto fespace = problem.getBlockAssembler()->block(mtn.second[0])->getFEspace(fieldIndex);

                std::vector<shp<VECTOREPETRA>> snapshots;

                while(std::getline(infile,line))
                {
                    shp<VECTOREPETRA> newVector(new VECTOREPETRA(fespace->map()));

                    std::stringstream linestream(line);
                    std::string value;
                    unsigned int j = 0;
                    while(getline(linestream,value,','))
                    {
                        newVector->operator[](j) = std::atof(value.c_str());
                        j++;
                    }
                    if (j != newVector->epetraVector().GlobalLength())
                        throw new Exception("Stored snapshot length does not match fespace dimension!");

                    snapshots.push_back(newVector);
                }

                unsigned int nsnapshots = snapshots.size() / mtn.second.size();

                for (unsigned int k = 0; k < nsnapshots; k++)
                {
                    unsigned int count = 0;
                    for (auto bindex : mtn.second)
                    {
                        shp<DistributedVector> vecWrap(new DistributedVector());
                        vecWrap->setVector(snapshots[k + count * nsnapshots]);
                        convert<BlockVector>(convert<BlockVector>(auxVec)->block(bindex))->setBlock(fieldIndex,vecWrap);
                        count++;
                    }

                    const bool inverse = true;
                    std::cout << "Applying piola" << std::endl << std::flush;
                    problem.getBlockAssembler()->applyPiola(auxVec, inverse);
                }

                infile.close();

                std::ofstream outfile;
                outfile.open(snapshotFile, omode);
                for (auto sol : snapshots)
                {
                    shp<DistributedVector> vectorWrap(new DistributedVector());
                    vectorWrap->setVector(sol);
                    std::string str2write = vectorWrap->getString(',') + "\n";
                    if (M_comm->MyPID() == 0)
                        outfile.write(str2write.c_str(), str2write.size());
                }
                outfile.close();
            }

        }
    }
}

std::vector<double>
SnapshotsSampler::
sampleParametersInflow()
{
    std::vector<double> array_params;

    std::string inflow_type = M_data("rb/offline/snapshots/inflow_type", "default");
    bool withOutflow = M_data("rb/offline/snapshots/add_outflow_param", false);
    int numInletConditions = M_data("bc_conditions/numinletbcs", 1);
    unsigned int numOutletConditions = M_data("bc_conditions/numoutletbcs", 0);

    std::vector<std::array<double,2>> param_bounds;

    if ((!std::strcmp(inflow_type.c_str(), "default")) || (!std::strcmp(inflow_type.c_str(), "periodic")))
    {
        param_bounds.push_back(std::array<double,2>({M_data("rb/offline/snapshots/a_min", 0.0),
                                                       M_data("rb/offline/snapshots/a_max", 1.0)}));
        param_bounds.push_back(std::array<double,2>({M_data("rb/offline/snapshots/c_min", 0.0),
                                                       M_data("rb/offline/snapshots/c_max", 1.0)}));
    }
    else if ((!std::strcmp(inflow_type.c_str(), "systolic")) || (!std::strcmp(inflow_type.c_str(), "heartbeat")))
    {
        param_bounds.push_back(std::array<double,2>({M_data("rb/offline/snapshots/Dt_min", 0.0),
                                                     M_data("rb/offline/snapshots/Dt_max", 1.0)}));
        param_bounds.push_back(std::array<double,2>({M_data("rb/offline/snapshots/DV0_min", 0.0),
                                                     M_data("rb/offline/snapshots/DV0_max", 1.0)}));
        param_bounds.push_back(std::array<double,2>({M_data("rb/offline/snapshots/DVM_min", 0.0),
                                                     M_data("rb/offline/snapshots/DVM_max", 1.0)}));
        param_bounds.push_back(std::array<double,2>({M_data("rb/offline/snapshots/DVm_min", 0.0),
                                                     M_data("rb/offline/snapshots/DVm_max", 1.0)}));
        if (!std::strcmp(inflow_type.c_str(), "heartbeat"))
            param_bounds.push_back(std::array<double,2>({M_data("rb/offline/snapshots/DVMd_min", 0.0),
                                                         M_data("rb/offline/snapshots/DVMd_max", 1.0)}));
    }

    unsigned int num_params_inflow = param_bounds.size();

    for (unsigned int numInlet=1; numInlet < numInletConditions; numInlet++)
    {
        std::string dataEntryMin = "rb/offline/snapshots/in_min_" + std::to_string(numInlet);
        std::string dataEntryMax = "rb/offline/snapshots/in_max_" + std::to_string(numInlet);
        param_bounds.push_back(std::array<double,2>({M_data(dataEntryMin, 0.0),
                                                     M_data(dataEntryMax, 1.0)}));
    }

    if (withOutflow)
        for (unsigned int numOutlet=0; numOutlet < numOutletConditions; numOutlet++)
        {
            std::string dataEntry = "bc_conditions/outlet" + std::to_string(numOutlet);
            std::string dataEntryMin = "rb/offline/snapshots/out_min_" + std::to_string(numOutlet);
            std::string dataEntryMax = "rb/offline/snapshots/out_max_" + std::to_string(numOutlet);
            if (!std::strcmp(M_data(dataEntry + "/type", "windkessel").c_str(), "dirichlet"))
                param_bounds.push_back(std::array<double,2>({M_data(dataEntryMin, 0.0),
                                                             M_data(dataEntryMax, 1.0)}));
        }

    std::vector<double> vec = sampleParameters(param_bounds);

    auto inletBC = std::bind(M_inflow,
                             std::placeholders::_1, vec);

    if (numInletConditions > 1)
    {
        double sum = 1.0 + std::accumulate(vec.begin()+num_params_inflow,
                                           vec.begin()+num_params_inflow+numInletConditions-1,
                                           0.0);
        for (auto it = vec.begin()+num_params_inflow; it != vec.begin()+num_params_inflow+numInletConditions-1; ++it)
            *it /= sum;
        vec.insert(vec.begin() + num_params_inflow, 1.0 / sum);
    }
    else
        vec.insert(vec.begin() + num_params_inflow, 1.0);

    array_params = vec;  // setting parameter values for exporting

    for (unsigned int numInlet=0; numInlet < numInletConditions; numInlet++)
    {
        unsigned int cnt = num_params_inflow + numInlet;
        std::function<double(double)> inletDirichlet = [vec, cnt, inletBC] (double t)
                {return vec[cnt] * inletBC(t);};
        M_data.setInletBC(inletDirichlet, numInlet);
    }

    if (withOutflow)
    {
        unsigned int cnt = num_params_inflow + numInletConditions;

        for (unsigned int numOutlet=0; numOutlet < numOutletConditions; numOutlet++)
        {
            std::string dataEntry = "bc_conditions/outlet" + std::to_string(numOutlet);
            if ((!std::strcmp(M_data(dataEntry + "/type", "windkessel").c_str(), "neumann")) ||
            (!std::strcmp(M_data(dataEntry + "/type", "windkessel").c_str(), "windkessel")) ||
            (!std::strcmp(M_data(dataEntry + "/type", "windkessel").c_str(), "coronary")))
                throw new Exception("Invalid outlet BC type! Only Dirichlet BCs are supported.");
            else if (!std::strcmp(M_data(dataEntry + "/type", "windkessel").c_str(), "dirichlet"))
            {
                std::function<double(double)> outletDirichlet = [vec, cnt, inletBC] (double t)
                        {return vec[cnt] * inletBC(t);};
                M_data.setOutletBC(outletDirichlet, numOutlet);
            }
            cnt++;
        }
    }

    return array_params;
}

std::vector<double>
SnapshotsSampler::
sampleParametersPhysics()
{
    std::vector<double> array_params;

    bool sampleFluidPhysics = M_data("rb/offline/snapshots/sample_fluid_physics", false);
    bool sampleStructurePhysics = M_data("rb/offline/snapshots/sample_structure_physics", false);
    bool sampleClothPhysics = M_data("rb/offline/snapshots/sample_cloth_physics", false);

    if (sampleFluidPhysics)
    {
        std::vector<std::array<double,2>> param_bounds;
        std::vector<double> param_ref_values;

        param_bounds.push_back(std::array<double,2>({M_data("rb/offline/snapshots/rho_f_min", -0.1),
                                                       M_data("rb/offline/snapshots/rho_f_max", 0.1)}));
        param_ref_values.push_back(1.06);
        param_bounds.push_back(std::array<double,2>({M_data("rb/offline/snapshots/mu_f_min", -0.1),
                                                       M_data("rb/offline/snapshots/mu_f_max", 0.5)}));
        param_ref_values.push_back(0.035);

        std::vector<double> array_params_fluid = sampleParameters(param_bounds);
        array_params.insert(array_params.end(), array_params_fluid.begin(), array_params_fluid.end());

        M_data.setValueDouble("fluid/density", (1.0+array_params_fluid[0]) * param_ref_values[0]);
        M_data.setValueDouble("fluid/viscosity", (1.0+array_params_fluid[1]) * param_ref_values[1]);
    }

    if (sampleStructurePhysics)
    {
        std::vector<std::array<double,2>> param_bounds;
        std::vector<double> param_ref_values;

        param_bounds.push_back(std::array<double,2>({M_data("rb/offline/snapshots/rho_s_min", 0.0),
                                                     M_data("rb/offline/snapshots/rho_s_max", 0.2)}));
        param_ref_values.push_back(1.2);
        param_bounds.push_back(std::array<double,2>({M_data("rb/offline/snapshots/h_s_min", -0.5),
                                                     M_data("rb/offline/snapshots/h_s_max", 1.0)}));
        param_ref_values.push_back(0.1);
        param_bounds.push_back(std::array<double,2>({M_data("rb/offline/snapshots/E_min", -0.5),
                                                     M_data("rb/offline/snapshots/E_max", 0.5)}));
        param_ref_values.push_back(4e6);
        param_bounds.push_back(std::array<double,2>({M_data("rb/offline/snapshots/nu_min", -0.4),
                                                     M_data("rb/offline/snapshots/nu_max", 0.0)}));
        param_ref_values.push_back(0.5);

        std::vector<double> array_params_structure = sampleParameters(param_bounds);
        array_params.insert(array_params.end(), array_params_structure.begin(), array_params_structure.end());

        M_data.setValueDouble("structure/density", (1.0+array_params_structure[0]) * param_ref_values[0]);
        M_data.setValueDouble("structure/thickness", (1.0+array_params_structure[1]) * param_ref_values[1]);
        M_data.setValueDouble("structure/constant_thickness", 1);  // setting constant thickness here !!
        M_data.setValueDouble("structure/young", (1.0+array_params_structure[2]) * param_ref_values[2]);
        M_data.setValueDouble("structure/poisson", (1.0+array_params_structure[3]) * param_ref_values[3]);
    }

    if (sampleClothPhysics)
    {
        std::vector<std::array<double,2>> param_bounds;

        unsigned int n_cloths = M_data("cloth/n_cloths", 0.0);
        for (unsigned int i=0; i<n_cloths; i++)
            param_bounds.push_back(std::array<double,2>({M_data("rb/offline/snapshots/Rcloth_min", 1.0),
                                                          M_data("rb/offline/snapshots/Rcloth_max", 3.0)}));

        std::vector<double> array_exp_cloth = sampleParameters(param_bounds, true);
        std::vector<double> array_params_cloth;
        array_params_cloth.resize(array_exp_cloth.size());
        std::transform(array_exp_cloth.begin(), array_exp_cloth.end(), array_params_cloth.begin(),
                       [](double x){return std::pow(10.0, x) * (x > 1e-8);});
        array_params.insert(array_params.end(), array_params_cloth.begin(), array_params_cloth.end());

        for (unsigned int i=0; i<n_cloths; i++)
        {
            std::string path = "cloth/cloth" + std::to_string(i) + "/density";
            M_data.setValueDouble(path, array_params_cloth[i]);
        }
    }

    return array_params;
}

std::vector<double>
SnapshotsSampler::
sampleParameters(const std::vector<std::array<double,2>>& param_bounds,
                 const bool is_cloth) const
{
    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<double> out;

    for (auto elem : param_bounds)
    {
        bool is_active;
        if (is_cloth)
        {
            unsigned int n_cloths = param_bounds.size();
            std::bernoulli_distribution bernoulli_distribution(1.0 / n_cloths);
            is_active = bernoulli_distribution(gen);
        }

        if ((!is_cloth) or (is_active))
        {
            std::uniform_real_distribution<> distribution(elem[0], elem[1]);
            out.push_back(distribution(gen));
        }
        else
            out.push_back(0.0);
    }

    return out;
}

}  // namespace RedMA

