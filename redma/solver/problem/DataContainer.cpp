#include "DataContainer.hpp"

namespace RedMA
{

DataContainer::
DataContainer()
{

}

void
DataContainer::
setDatafile(const std::string& datafile)
{
    M_datafile.reset(new GetPot(datafile));
}

void
DataContainer::
setInflow(const std::function<double(double)>& inflow)
{
    M_inflow = inflow;
}

void
DataContainer::
setInflowDt(const std::function<double(double)>& inflowDt)
{
    M_inflowDt = inflowDt;
}

void
DataContainer::
setVerbose(bool verbose)
{
    M_verbose = verbose;

}

std::string
DataContainer::
operator()(std::string location, std::string defValue) const
{
    return M_datafile->operator()(location.c_str(), defValue.c_str());
}

int
DataContainer::
operator()(std::string location, int defValue) const
{
    return M_datafile->operator()(location.c_str(), defValue);
}

double
DataContainer::
operator()(std::string location, double defValue) const
{
    return M_datafile->operator()(location.c_str(), defValue);
}

}
