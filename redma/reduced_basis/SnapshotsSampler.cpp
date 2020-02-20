#include "SnapshotsSampler.hpp"

namespace RedMA
{

SnapshotsSampler::
SnapshotsSampler(const DataContainer& data, EPETRACOMM comm) :
  M_problem(data, comm, false),
  M_data(data),
  M_comm(comm)
{
}

void
SnapshotsSampler::
takeSnapshots()
{
    boost::filesystem::create_directory("snapshots");

    GeometryPrinter printer;

    unsigned int nSnapshots = M_data("snapshots/number", 10);
    unsigned int bound = M_data("snapshots/bound", 0.2);

    for (unsigned int i = 0; i < nSnapshots; i++)
    {
        M_problem.getTree().randomSampleAroundOriginalValue(bound);
        M_problem.setup();
        M_problem.solve();

        std::string filename = "snapshots/snapshots_tree" + std::to_string(i) + ".xml";
        printer.saveToFile(M_problem.getTree(), filename, M_comm);
    }
}

}  // namespace RedMA
