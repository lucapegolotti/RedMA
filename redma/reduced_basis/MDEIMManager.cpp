#include "MDEIMManager.hpp"

namespace RedMA
{

MDEIMManager::
MDEIMManager(const DataContainer& dataContainer, EPETRACOMM comm,
             std::map<unsigned int, std::string> idmeshmap) :
  M_data(dataContainer),
  M_comm(comm),
  M_IDMeshMap(idmeshmap)
{
}

void
MDEIMManager::
load()
{
    Chrono chrono;
    chrono.start();

    std::string msg = "[MDEIMManager] loading structures ...";
    printlog(YELLOW, msg, this->M_data.getVerbose());

    using namespace boost::filesystem;

    std::string mdeimdir = M_data("rb/online/mdeim/directory", "mdeims");

    for (auto idmeshpair : M_IDMeshMap)
    {
        std::string curdir = mdeimdir + "/" + idmeshpair.second;
        std::vector<SHP(BlockMDEIM)> mdeimvec;

        unsigned int matrixIndex = 0;
        while (exists(curdir + "/blockmdeim" + std::to_string(matrixIndex)))
        {
            std::string blockmdeimdir = curdir + "/blockmdeim" + std::to_string(matrixIndex);

            unsigned int countrow = 0;
            while (exists(blockmdeimdir + "/mdeim_" + std::to_string(countrow) + "_0"))
                countrow++;

            unsigned int countcol = 0;
            while (exists(blockmdeimdir + "/mdeim_0_" + std::to_string(countcol)))
                countcol++;

            SHP(BlockMDEIM) newMdeim(new BlockMDEIM());
            newMdeim->setDataContainer(M_data);

            newMdeim->setComm(M_comm);
            newMdeim->resize(countrow, countcol);
            newMdeim->setMatrixIndex(matrixIndex);
            newMdeim->loadMDEIM(blockmdeimdir);
            mdeimvec.push_back(newMdeim);
            matrixIndex++;
        }
        M_mdeims[idmeshpair.second] = mdeimvec;
    }

    msg = "done, in ";
    msg += std::to_string(chrono.diff());
    msg += " seconds\n";
    printlog(YELLOW, msg, this->M_data.getVerbose());
}

}  // namespace RedMA
