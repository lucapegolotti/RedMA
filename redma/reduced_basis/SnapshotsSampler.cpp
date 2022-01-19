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
    bool withOutflow = M_data("rb/offline/snapshots/add_outflow_param", false);
    int numInletConditions = M_data("bc_conditions/numinletbcs", -1);
    unsigned int numOutletConditions = M_data("bc_conditions/numoutletbcs", 0);
    if (numInletConditions == -1)
        numInletConditions = 1;

    fs::create_directory(outdir);
    GeometryPrinter printer;

    unsigned int nSnapshots = M_data("rb/offline/snapshots/number", 10);
    double bound = M_data("rb/offline/snapshots/bound", 0.2);
    std::vector<double> array_params;

    for (unsigned int i = 0; i < nSnapshots; i++)
    {
        // to guarantee (almost...) that two snapshots are not saved at the same location!
        unsigned int paramIndex = Nstart;
        
        // we find the first parameter index available, starting from Nstart
        while (fs::exists(outdir + "/param" + std::to_string(paramIndex)))
            paramIndex++;
        std::string curdir = outdir + "/param" + std::to_string(paramIndex);

        if (!std::strcmp(param_type.c_str(), "inflow"))
        {
            std::vector<std::vector<double>> param_bounds;

            for (unsigned int numInlets=0; numInlets < numInletConditions; numInlets++)
            {
                param_bounds.push_back(std::vector<double>({M_data("rb/offline/snapshots/a_min", 0.0),
                                                            M_data("rb/offline/snapshots/a_max", 1.0)}));
                param_bounds.push_back(std::vector<double>({M_data("rb/offline/snapshots/c_min", 0.0),
                                                            M_data("rb/offline/snapshots/c_max", 1.0)}));
            }

            if (withOutflow)
                for (unsigned int numOutlet=0; numOutlet < numOutletConditions; numOutlet++)
                {
                    std::string dataEntry = "bc_conditions/outlet" + std::to_string(numOutlet);
                    if (!std::strcmp(M_data(dataEntry + "/type", "windkessel").c_str(), "dirichlet"))
                        param_bounds.push_back(std::vector<double>({M_data("rb/offline/snapshots/out_min", 0.0),
                                                                      M_data("rb/offline/snapshots/out_max", 1.0)}));
                }

            std::vector<double> vec = inflowSnapshots(param_bounds);
            array_params = vec;

            for (unsigned int numInlets=0; numInlets < numInletConditions; numInlets++)
            {
                auto inletBC = std::bind(M_inflow,
                                         std::placeholders::_1,
                                         vec[2*numInlets],vec[2*numInlets+1]);
                M_data.setInletBC(inletBC, numInlets);
            }

            if (withOutflow)
            {
                unsigned int cnt = 2*numInletConditions;
                auto inletBC = std::bind(M_inflow,
                                         std::placeholders::_1,
                                         vec[0],vec[1]);
                for (unsigned int numOutlet = 0; numOutlet < numOutletConditions; numOutlet++)
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
                }
            }
        }

        GlobalProblem problem(M_data, M_comm, false);

        problem.doStoreSolutions();

        fs::create_directory(curdir);

        if (!std::strcmp(param_type.c_str(), "geometric"))
            problem.getTree().randomSampleAroundOriginalValue(bound);

        problem.setup();

        if (!problem.isFEProblem())
            throw new Exception("The tree must be composed of only FE nodes to "
                                "sample the snapshots!");

        std::string filename = curdir + "/tree.xml";
        printer.saveToFile(problem.getTree(), filename, M_comm);
        
	clock_t t1,t2;
	t1 = clock();

        problem.solve();

	t2 = clock();
	double time_dif = (double)(t2 - t1)/CLOCKS_PER_SEC;	

        dumpSnapshots(problem, curdir, time_dif, array_params);
    }

}

void
SnapshotsSampler::
dumpSnapshots(GlobalProblem& problem,
              std::string outdir,
              const double computational_time,
              const std::vector<double> array_params = {})
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
    double density = M_data("rb/offline/fluid/density", 1.0);
    double viscosity = M_data("rb/offline/fluid/viscosity", 1.0);

    std::ios_base::openmode omode = std::ios_base::app;
    if (binary)
        omode = omode | std::ios::binary;

    M_mass->block(0,0)->dump("M");
    M_stiffness->block(0,0)->dump("A");
    M_divergence->block(0,1)->dump("BdivT");
    M_divergence->block(1,0)->dump("Bdiv");

    /*for (auto sol : solutions)
        problem.getBlockAssembler()->applyPiola(sol, true);*/

    for (auto idmeshtype : IDmeshTypeMap)
    {
        std::string meshtypedir = outdir + "/" + idmeshtype.second;
        fs::create_directory(meshtypedir);

        unsigned int nfields = solutions[0]->block(idmeshtype.first)->nRows();

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


        for(unsigned int j = 0; j < n_dual_blocks; j++)  // save Lagrange multipliers, if any

        {
            std::string outfilename = meshtypedir + "/lagmult_" + std::to_string(j) + ".snap";
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

        std::ofstream file(outdir + "/computational_time.txt", std::ios_base::app);
	file << std::fixed << std::setprecision(2) << computational_time << std::endl;

    }
}

void
SnapshotsSampler::
transformSnapshotsWithPiola(std::string snapshotsDir,
                            unsigned int fieldIndex,
                            unsigned int maxSnapshot)
{

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
inflowSnapshots(const std::vector<std::vector<double>>& param_bounds)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<double> out;

    for (auto elem : param_bounds)
    {
        std::uniform_real_distribution<> distribution(elem[0], elem[1]);
        out.push_back(distribution(gen));
    }

    return out;
}

}  // namespace RedMA

