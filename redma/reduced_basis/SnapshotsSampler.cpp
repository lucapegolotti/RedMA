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
takeSnapshots()
{
    // using namespace boost::filesystem;
    std::string outdir = M_data("rb/offline/snapshots/directory", "snapshots");
    // if (fs::exists(outdir))
    //    throw new Exception("Snapshots directory already exists!");
    fs::create_directory(outdir);
    GeometryPrinter printer;

    unsigned int nSnapshots = M_data("rb/offline/snapshots/number", 10);
    double bound = M_data("rb/offline/snapshots/bound", 0.2);

    for (unsigned int i = 0; i < nSnapshots; i++)
    {
        GlobalProblem problem(M_data, M_comm, false);
        // this is to allow taking the snapshots at the end of the simulation from
        // the fem problem
        problem.doStoreSolutions();

        unsigned int paramIndex = 0;
        while (fs::exists(outdir + "/param" + std::to_string(paramIndex)))
            paramIndex++;
        std::string curdir = outdir + "/param" + std::to_string(paramIndex);
        fs::create_directory(curdir);
        problem.getTree().randomSampleAroundOriginalValue(bound);
        problem.setup();

        std::string filename = curdir + "/tree.xml";
        printer.saveToFile(problem.getTree(), filename, M_comm);

        problem.solve();
        dumpSnapshots(problem, curdir);
    }
}

void
SnapshotsSampler::
dumpSnapshots(GlobalProblem& problem,
              std::string outdir)
{
    auto IDmeshTypeMap = problem.getBlockAssembler()->getIDMeshTypeMap();
    auto solutions = problem.getSolutions();

    unsigned int takeEvery = M_data("rb/offline/snapshots/take_every", 5);
    bool binary = M_data("rb/offline/snapshots/dumpbinary", true);
    bool computereynolds = M_data("rb/offline/snapshots/computereynolds", true);
    double density = M_data("rb/offline/fluid/density", 1.0);
    double viscosity = M_data("rb/offline/fluid/viscosity", 1.0);

    std::ios_base::openmode omode = std::ios_base::app;
    if (binary)
        omode = omode | std::ios::binary;

    for (auto sol : solutions)
        problem.getBlockAssembler()->applyGlobalPiola(sol, true);

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
                auto solBlck = convert<DistributedVector>(convert<BlockVector>(convert<BlockVector>(sol)->block(idmeshtype.first))->block(i));
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
            reynoldsfile.open(meshtypedir + "/reynolds.txt", std::ios_base::app |
                                                             std::ios::binary);
            for (auto sol : solutions)
            {
                auto solBlck = convert<DistributedVector>(convert<BlockVector>(convert<BlockVector>(sol)->block(idmeshtype.first))->block(0));
                double Umax = solBlck->maxMagnitude3D();
                auto tNode = problem.getBlockAssembler()->block(0)->getTreeNode();
                double D = 2 * tNode->M_block->getInlet().M_radius;
                double curReynolds = Umax * density * D / viscosity;

                if (M_comm->MyPID() == 0)
                    reynoldsfile << std::to_string(curReynolds) << std::endl;
            }
            reynoldsfile.close();
        }
    }
}

void
SnapshotsSampler::
transformSnapshotsWithPiola(std::string snapshotsDir,
                            unsigned int fieldIndex,
                            unsigned int maxSnapshot)
{
    // using namespace boost::filesystem;

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

            // for each mesh type, assemblers sharing that mesh
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
                    unsigned int i = 0;
                    while(getline(linestream,value,','))
                    {
                        newVector->operator[](i) = std::atof(value.c_str());
                        i++;
                    }
                    if (i != newVector->epetraVector().GlobalLength())
                        throw new Exception("Stored snapshot length does not match fespace dimension!");

                    snapshots.push_back(newVector);
                }

                unsigned int nsnapshots = snapshots.size() / mtn.second.size();

                for (unsigned int i = 0; i < nsnapshots; i++)
                {
                    unsigned int count = 0;
                    for (auto bindex : mtn.second)
                    {
                        shp<DistributedVector> vecWrap(new DistributedVector());
                        vecWrap->setVector(snapshots[i + count * nsnapshots]);
                        convert<BlockVector>(convert<BlockVector>(auxVec)->block(bindex))->setBlock(fieldIndex,vecWrap);
                        count++;
                    }

                    const bool inverse = true;
                    std::cout << "Applying piola" << std::endl << std::flush;
                    problem.getBlockAssembler()->applyGlobalPiola(auxVec, inverse);
                }

                infile.close();

                std::ofstream outfile;
                outfile.open(snapshotFile, omode);
                unsigned int count = 0;
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

}  // namespace RedMA
