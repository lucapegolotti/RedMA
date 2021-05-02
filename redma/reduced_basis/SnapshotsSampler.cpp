#include "SnapshotsSampler.hpp"

namespace RedMA
{

SnapshotsSampler::
SnapshotsSampler(const DataContainer& data, EPETRACOMM comm, bool geo) :
  M_data(data),
  M_comm(comm),
  geometry(geo)
{
}


SnapshotsSampler::
SnapshotsSampler(const DataContainer& data, EPETRACOMM comm, bool geo ,myFunctionType inflow_1, myFunctionType inflow_2,
                 double index_,double index_max_) :
    M_data(data),
    M_comm(comm),
    geometry(geo),
    analyticInflow(true),
    inflow1(inflow_1),inflow2(inflow_2),
    index(index_),index_max(index_max_)
    {
    }

void
SnapshotsSampler::
takeSnapshots()
{
    std::string outdir = M_data("rb/offline/snapshots/directory", "snapshots");

    fs::create_directory(outdir);
    GeometryPrinter printer;

    unsigned int nSnapshots = M_data("rb/offline/snapshots/number", 10);

    double bound = M_data("rb/offline/snapshots/bound", 0.2);
    if(geometry){
    for (unsigned int i = 0; i < nSnapshots; i++)
    {


            unsigned int paramIndex = 0;
        // we find the first parameter index available
            while (fs::exists(outdir + "/param" + std::to_string(paramIndex)))
                paramIndex++;

            std::string curdir = outdir + "/param" + std::to_string(paramIndex);
            GlobalProblem problem(M_data, M_comm, false);
            // this is to allow taking the snapshots at the end of the simulation from
            // the fem problem
            problem.doStoreSolutions();

            fs::create_directory(curdir);
            problem.getTree().randomSampleAroundOriginalValue(bound);
            problem.setup();

            if (!problem.isFEProblem())
                throw new Exception("The tree must be composed of only FE nodes to "
                                    "sample the snapshots!");

            std::string filename = curdir + "/tree.xml";
            printer.saveToFile(problem.getTree(), filename, M_comm);

            problem.solve();
            dumpSnapshots(problem, curdir);}

        }else {


        double alpha_min=0;
        double alpha_max=1;
        double index_min=1;
        double alpha=alpha_min+(index-index_min-1)/(index_max-index_min-1)*(alpha_max-alpha_min);
        if (analyticInflow) {
            auto inflow1_alpha([alpha, this](double t) {
                return alpha * inflow1(t);
             });
             auto inflow2_alpha([alpha, this](double t) {
                return (1 - alpha) * inflow2(t);
            });
            M_data.setInflow(inflow1_alpha, 2);
            M_data.setInflow(inflow2_alpha, 3);
        } else {

        }
         std::string curdir = outdir + "/param" + std::to_string(alpha);
        GlobalProblem problem(M_data, M_comm, false);
        // this is to allow taking the snapshots at the end of the simulation from
        // the fem problem
        problem.doStoreSolutions();

        fs::create_directory(curdir);
        // problem.getTree().randomSampleAroundOriginalValue(bound);
        problem.setup();

        if (!problem.isFEProblem())
            throw new Exception("The tree must be composed of only FE nodes to "
                                    "sample the snapshots!");

        std::string filename = curdir + "/tree.xml";
        printer.saveToFile(problem.getTree(), filename, M_comm);

        problem.solve();
        dumpSnapshots(problem, curdir);}
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

    if (geometry)
    {
        for (auto sol : solutions)
            problem.getBlockAssembler()->applyPiola(sol, true);
    }

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
                double D = 2 * tNode->M_block->getInlet(0).M_radius;
                double curReynolds = Umax * density * D / viscosity;

                if (M_comm->MyPID() == 0)
                    reynoldsfile << std::to_string(curReynolds) << std::endl;
            }
            reynoldsfile.close();
        }
    }
}

}  // namespace RedMA
