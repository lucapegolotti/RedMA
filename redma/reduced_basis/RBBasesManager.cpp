#include "RBBasesManager.hpp"

namespace RedMA
{

RBBasesManager::
RBBasesManager(const DataContainer& dataContainer, EPETRACOMM comm,
             std::map<unsigned int, std::string> idmeshmap) :
  M_data(dataContainer),
  M_comm(comm),
  M_IDMeshMap(idmeshmap)
{
}

void
RBBasesManager::
load()
{
    LifeV::LifeChrono chrono;
    chrono.start();

    std::string msg = "[RBBasesManager] loading singular values ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    using namespace boost::filesystem;

    std::string basisdir = M_data("rb/online/basis/directory", "basis");

    for (auto idmeshpair : M_IDMeshMap)
    {
        std::string curdir = basisdir + "/" + idmeshpair.second;

        // count number of fields

        unsigned int numFields = 0;
        while (exists(curdir + "/svd" + std::to_string(numFields) + ".txt"))
            numFields++;

        SHP(RBBases) newBases(new RBBases(M_data, M_comm));
        newBases->setNumberOfFields(numFields);
        newBases->setPath(curdir);
        newBases->loadSingularValues();
        M_bases[idmeshpair.second] = newBases;
    }

    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());
}

void
RBBasesManager::
loadBases()
{
    LifeV::LifeChrono chrono;
    chrono.start();

    std::string msg = "[RBBasesManager] loading bases ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    for (auto& bases : M_bases)
        bases.second->loadBases();

    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());
}

}  // namespace RedMA
