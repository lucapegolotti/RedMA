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
loadSingularValues()
{
    Chrono chrono;
    chrono.start();

    std::string msg = "[RBBasesManager] loading singular values ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    // using namespace boost::filesystem;

    std::string basisdir = M_data("rb/online/basis/directory", "basis");

    for (auto idmeshpair : M_IDMeshMap)
    {
        std::string curdir = basisdir + "/" + idmeshpair.second;

        if (fs::exists(curdir))
        {
            // count number of fields
            unsigned int numFields = 0;
            while (fs::exists(curdir + "/svd" + std::to_string(numFields) + ".txt"))
                numFields++;

            shp<RBBases> newBases(new RBBases(M_data, M_comm));
            newBases->setNumberOfFields(numFields);
            newBases->setPath(curdir);
            newBases->loadSingularValues();
            M_bases[idmeshpair.second] = newBases;
        }
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
    Chrono chrono;
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

shp<RBBases>
RBBasesManager::
getRBBases(std::string meshName)
{
    if (M_bases.find(meshName) == M_bases.end())
        throw new Exception("Basis with mesh " + meshName + "is not available!");

    return M_bases[meshName];
}

}  // namespace RedMA
