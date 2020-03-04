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
    using namespace boost::filesystem;

    std::string outdir = M_data("rb/offline/snapshots/directory", "snapshots");

    if (exists(outdir))
        throw new Exception("Snapshots directory already exists!");

    create_directory(outdir);

    GeometryPrinter printer;

    unsigned int nSnapshots = M_data("rb/offline/snapshots/number", 10);
    double bound = M_data("rb/offline/snapshots/bound", 0.2);

    for (unsigned int i = 0; i < nSnapshots; i++)
    {
        ProblemFEM problem(M_data, M_comm, false);
        // this is to allow taking the snapshots at the end of the simulation from
        // the fem problem
        problem.doStoreSolutions();

        std::string curdir = outdir + "/param" + std::to_string(i);
        create_directory(curdir);
        problem.getTree().randomSampleAroundOriginalValue(bound);
        problem.setup();
        problem.solve();

        std::string filename = curdir + "/tree.xml";
        printer.saveToFile(problem.getTree(), filename, M_comm);

        dumpSnapshots(problem, curdir);
    }
}

void
SnapshotsSampler::
dumpSnapshots(ProblemFEM& problem,
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

    for (auto idmeshtype : IDmeshTypeMap)
    {
        std::string meshtypedir = outdir + "/" + idmeshtype.second;
        boost::filesystem::create_directory(meshtypedir);

        unsigned int nfields = solutions[0].block(idmeshtype.first).nRows();

        for (unsigned int i = 0; i < nfields; i++)
        {
            std::string outfilename = meshtypedir + "/field" + std::to_string(i) + ".snap";

            std::ofstream outfile;
            outfile.open(outfilename, omode);
            unsigned int count = 0;
            for (auto sol : solutions)
            {
                if (count % takeEvery == 0)
                {
                    std::string str2write = sol.block(idmeshtype.first).block(i).getString(',') + "\n";
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
                double Umax = sol.block(idmeshtype.first).block(0).maxMagnitude3D();
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

}  // namespace RedMA
