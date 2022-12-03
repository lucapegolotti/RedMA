#include <redma/reduced_basis/aMatricesGenerator.hpp>

namespace RedMA
{

aMatricesGenerator::
aMatricesGenerator(const DataContainer& data, EPETRACOMM comm) :
M_data(data),
M_comm(comm)
{
    if (M_comm->MyPID() != 0)
        throw new Exception("MatricesGenerator does not support more than one proc");
}

void
aMatricesGenerator::
setDefaultParameterValues(const std::map<std::string, bool> &categories)
{
    auto it_fluid = categories.find("fluid");
    if ((it_fluid != categories.end()) && (it_fluid->second))
    {
        M_data.setValueDouble("fluid/density", 1.06);
        M_data.setValueDouble("fluid/viscosity", 0.035);
    }

    auto it_structure = categories.find("structure");
    if ((it_structure != categories.end()) && (it_structure->second))
    {
        M_data.setValueDouble("structure/density", 1.2);
        M_data.setValueDouble("structure/thickness", 0.1);
        M_data.setValueInt("structure/constant_thickness", 1);  // enforce constant thickness in RB
        M_data.setValueDouble("structure/young", 4e6);
        M_data.setValueDouble("structure/poisson", 0.5);
    }

    auto it_cloth = categories.find("cloth");
    if ((it_cloth != categories.end()) && (it_cloth->second))
    {
        unsigned int n_cloths = M_data("cloth/n_cloths", 0);
        for (unsigned int i=0; i<n_cloths; i++)
            M_data.setValueDouble("cloth/cloth" + std::to_string(i) + "/density", 1.0);
    }
}

void
aMatricesGenerator::
setDummyFlows()
{
    bool withOutflow = M_data("rb/offline/snapshots/add_outflow_param", false);
    int numInletConditions = M_data("bc_conditions/numinletbcs", 1);
    unsigned int numOutletConditions = M_data("bc_conditions/numoutletbcs", 0);

    std::function<double(double)> dummyFlow = [] (double t) {return 1.0;};

    for (unsigned int numInlet=0; numInlet < numInletConditions; numInlet++)
        M_data.setInletBC(dummyFlow, numInlet);

    if (withOutflow)
    {
        for (unsigned int numOutlet=0; numOutlet < numOutletConditions; numOutlet++)
        {
            std::string dataEntry = "bc_conditions/outlet" + std::to_string(numOutlet);
            if ((!std::strcmp(M_data(dataEntry + "/type", "windkessel").c_str(), "neumann")) ||
            (!std::strcmp(M_data(dataEntry + "/type", "windkessel").c_str(), "windkessel")) ||
            (!std::strcmp(M_data(dataEntry + "/type", "windkessel").c_str(), "coronary")))
                throw new Exception("Invalid outlet BC type! Only Dirichlet BCs are supported.");
            else if (!std::strcmp(M_data(dataEntry + "/type", "windkessel").c_str(), "dirichlet"))
                M_data.setOutletBC(dummyFlow, numOutlet);
        }
    }
}

}// namespace RedMA

