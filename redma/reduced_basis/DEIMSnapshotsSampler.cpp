#include "DEIMSnapshotsSampler.hpp"

namespace RedMA
{

DEIMSnapshotsSampler::
DEIMSnapshotsSampler(const DataContainer& data, EPETRACOMM comm) :
  M_data(data),
  M_comm(comm)
{
}

void
DEIMSnapshotsSampler::
takeDEIMSnapshots()
{
    using namespace boost::filesystem;

    std::string outdir = M_data("rb/offline/deim/directory", "deimsnapshots");
    create_directory(outdir);

    GeometryPrinter printer;

    unsigned int nSnapshots = M_data("rb/offline/deim/numbersnapshots", 10);
    double bound = M_data("rb/offline/deim/bound", 0.2);

    for (unsigned int i = 0; i < nSnapshots; i++)
    {
        ProblemRB problem(M_data, M_comm, false);
        problem.doStoreNonLinearTerms();
        problem.doStoreSolutions();

        unsigned int paramIndex = 0;
        while (exists(outdir + "/param" + std::to_string(paramIndex)))
            paramIndex++;
        std::string curdir = outdir + "/param" + std::to_string(paramIndex);
        create_directory(curdir);

        problem.getTree().randomSampleAroundOriginalValue(bound);
        problem.setup();
        problem.solve();

        std::string filename = curdir + "/tree.xml";
        printer.saveToFile(problem.getTree(), filename, M_comm);

        dumpDEIMSnapshots(problem, curdir);
    }
}

void
DEIMSnapshotsSampler::
dumpDEIMSnapshots(ProblemRB& problem,
                  std::string outdir)
{
    auto IDmeshTypeMap = problem.getBlockAssembler()->getIDMeshTypeMap();
    auto nonLinearTerms = problem.getNonLinearTerms();
    auto solutions = problem.getSolutions();

    unsigned int takeEvery = M_data("rb/offline/deim/take_snapshot_every", 5);
    bool binary = M_data("rb/offline/deim/dumpbinary_snapshots", true);

    std::ios_base::openmode omode = std::ios_base::app;
    if (binary)
        omode = omode | std::ios::binary;

    for (auto idmeshtype : IDmeshTypeMap)
    {
        std::string meshtypedir = outdir + "/" + idmeshtype.second;
        boost::filesystem::create_directory(meshtypedir);

        unsigned int nfields = nonLinearTerms[0].block(idmeshtype.first).nRows();

        for (unsigned int i = 0; i < nfields; i++)
        {
            std::string outfilename = meshtypedir + "/deimterm_fem" + std::to_string(i) + ".snap";

            std::ofstream outfiledeimfem;
            outfiledeimfem.open(outfilename, omode);

            outfilename = meshtypedir + "/deimterm_rb" + std::to_string(i) + ".snap";

            std::ofstream outfiledeimrb;
            outfiledeimrb.open(outfilename, omode);

            outfilename = meshtypedir + "/field_rb" + std::to_string(i) + ".snap";

            std::ofstream outfilerb;
            outfilerb.open(outfilename, omode);

            unsigned int count = 0;
            for (auto nlt : nonLinearTerms)
            {
                if (count % takeEvery == 0)
                {
                    auto nltFEM = problem.getBlockAssembler()->convertFunctionRBtoFEM(nlt, M_comm);
                    std::string str2write = nltFEM.block(idmeshtype.first).block(i).getString(',') + "\n";
                    outfiledeimfem.write(str2write.c_str(), str2write.size());

                    str2write = nlt.block(idmeshtype.first).block(i).getString(',') + "\n";
                    outfiledeimrb.write(str2write.c_str(), str2write.size());
                }
                count++;
            }

            count = 0;
            for (auto sol : solutions)
            {
                if (count % takeEvery == 0)
                {
                    std::string str2write = sol.block(idmeshtype.first).block(i).getString(',') + "\n";
                    outfilerb.write(str2write.c_str(), str2write.size());
                }
                count++;
            }
            outfiledeimfem.close();
            outfiledeimrb.close();
            outfilerb.close();
        }
    }
}

}  // namespace RedMA
