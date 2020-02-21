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
    std::string outdir = M_data("snapshots/directory", "snapshots");
    boost::filesystem::create_directory(outdir);

    GeometryPrinter printer;

    unsigned int nSnapshots = M_data("snapshots/number", 10);
    unsigned int bound = M_data("snapshots/bound", 0.2);

    for (unsigned int i = 0; i < nSnapshots; i++)
    {
        ProblemFEM problem(M_data, M_comm, false);
        // this is to allow taking the snapshots at the end of the simulation from
        // the fem problem
        problem.doStoreSolutions();

        std::string curdir = outdir + "/param" + std::to_string(i);
        boost::filesystem::create_directory(curdir);
        problem.getTree().randomSampleAroundOriginalValue(bound);
        problem.setup();
        problem.solve();

        std::string filename = curdir + "/tree.xml";
        printer.saveToFile(problem.getTree(), filename, M_comm);

        auto IDmeshTypeMap = problem.getBlockAssembler()->getIDMeshTypeMap();
        auto solutions = problem.getSolutions();
        dumpSnapshots(solutions, IDmeshTypeMap, curdir);
    }
}

void
SnapshotsSampler::
dumpSnapshots(const std::vector<BlockVector<BlockVector<FEVECTOR>>>& solutions,
              std::map<unsigned int, std::string> IDmeshTypeMap,
              std::string outdir)
{
    unsigned int takeEvery = M_data("snapshots/take_every", 5);
    bool binary = M_data("snapshots/dumpbinary", true);
    std::ios_base::openmode omode = std::ios_base::app;
    if (binary)
        omode = omode | std::ios::binary;

    for (auto idmeshtype : IDmeshTypeMap)
    {
        unsigned int dashpos = idmeshtype.second.find("/");
        unsigned int formatpos = idmeshtype.second.find(".mesh");
        std::string actualmeshname = idmeshtype.second.substr(dashpos + 1, formatpos - dashpos - 1);

        std::string meshtypedir = outdir + "/" + actualmeshname;
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
                    if (binary)
                        outfile.write(str2write.c_str(), str2write.size());
                    else
                        outfile << str2write;
                }

                count++;
            }
            outfile.close();
        }
    }
}

}  // namespace RedMA
