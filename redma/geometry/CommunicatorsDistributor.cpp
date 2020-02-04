#include "CommunicatorsDistributor.hpp"

namespace RedMA
{

CommunicatorsDistributor::
CommunicatorsDistributor(const GetPot& datafile, commPtr_Type comm) :
  M_datafile(datafile),
  M_comm(comm)
{
    createCommunicators();

    CoutRedirecter ct;
    ct.redirect();
    M_parser.reset(new GeometryParser(datafile,
                                      datafile("geometric_structure/xmlfile","tree.xml"),
                                      M_comm, false));

    M_tree = M_parser->getTree();

    std::string geometriesDir = datafile("geometric_structure/geometries_dir",
                                         "../../../meshes/");

    M_tree.readMeshes(geometriesDir);
    ct.restore();
    std::map<unsigned int, unsigned int> nnodes = M_tree.getNumPointsMeshes();

    for (auto nn : nnodes)
        std::cout << nn.first << " " << nn.second << std::endl << std::flush;

    loadBalancing(nnodes);
}

void
CommunicatorsDistributor::
loadBalancing(std::map<unsigned int, unsigned int> npointsmeshes)
{
    unsigned int numProc = M_comm->NumProc();

    // trasform npointsmeshes in pairs to order it
    std::vector<std::pair<unsigned int, unsigned int>> pairs;
    for (auto itr = npointsmeshes.begin(); itr != npointsmeshes.end(); ++itr)
        pairs.push_back(*itr);

    sort(pairs.begin(), pairs.end(), [=](std::pair<unsigned int, unsigned int>& a,
                                         std::pair<unsigned int, unsigned int>& b)
    {
        return a.second > b.second;
    }
    );

    std::vector<unsigned int> nnodesperproc(numProc,0);

    // greedy approach. We always give the new building block to the process that
    // has fewer nodes. We want to avoid building blocks being shared among
    // processes
    unsigned count = 0;
    for (auto np : pairs)
    {
        unsigned int indexmin = -1;
        unsigned int minnodes = -1; // wraps around as unsigned int
        unsigned int count = 0;
        for (auto curnodes : nnodesperproc)
        {
            if (curnodes < minnodes)
            {
                minnodes = curnodes;
                indexmin = count;
            }
            count++;
        }

        nnodesperproc[indexmin] += np.second;
        std::cout << nnodesperproc[indexmin] << std::endl << std::flush;
        M_processMap[np.first] = indexmin;
    }

    count = 0;
    for (auto curnodes : nnodesperproc)
    {
        std::cout << "id = " << count << " nodes = " << curnodes << std::endl;
        count++;
    }

    for (auto pmap : M_processMap)
    {
        std::cout << "block id = " << pmap.first << " proc id = " << pmap.second << std::endl << std::flush;

    }
}

void
CommunicatorsDistributor::
createCommunicators()
{
    unsigned int numProc = M_comm->NumProc();
    M_communicators.resize(numProc, nullptr);

    MPI_Comm comm = (std::dynamic_pointer_cast<Epetra_MpiComm>(M_comm))->Comm();

    // Group initialization
    MPI_Group commGroup;
    MPI_Comm_group(comm, &commGroup);

    for (unsigned int i = 0; i < numProc; i++)
    {
        int* ranks = new int[1];
        ranks[0] = i;
        // Create parallel group
        MPI_Group parallelGroup;
        MPI_Group_incl(commGroup, 1, ranks, &parallelGroup);

        // Create parallel comm
        MPI_Comm newComm;
        MPI_Comm_create(comm, parallelGroup, &newComm);
        if (i == M_comm->MyPID())
        {
            M_communicators[i].reset(new Epetra_MpiComm(newComm));
        }
        delete[] ranks;
    }
    M_comm->Barrier();
}

}
