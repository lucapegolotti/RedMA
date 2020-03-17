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

std::function<double(double)>
DataContainer::
getDistalPressure(const unsigned int& outletIndex) const
{
    auto it = M_distalPressures.find(outletIndex);

    if (it == M_distalPressures.end())
        throw new Exception("Error in DataContainer: requested Distal pressure not present");

    return it->second;
}

void
DataContainer::
setDistalPressure(const std::function<double(double)>& pressure,
                  const unsigned int& indexOutlet)
{
    M_distalPressures[indexOutlet] = pressure;
}

void
DataContainer::
setVerbose(bool verbose)
{
    M_verbose = verbose;

}

std::string
DataContainer::
operator()(std::string location, const char* defValue) const
{
    return M_datafile->operator()(location.c_str(), defValue);
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

bool
DataContainer::
operator()(std::string location, bool defValue) const
{
    return M_datafile->operator()(location.c_str(), defValue);
}

void
DataContainer::
setValue(std::string location, std::string value)
{
    M_datafile->set(location.c_str(), value.c_str());
}

void
DataContainer::
setValue(std::string location, int value)
{
    M_datafile->set(location.c_str(), value);
}

void
DataContainer::
setValue(std::string location, double value)
{
    M_datafile->set(location.c_str(), value);
}

void
DataContainer::
setValue(std::string location, bool value)
{
    M_datafile->set(location.c_str(), value);
}

void
DataContainer::
finalize()
{
    if (!M_inflow)
        generateInflow();
}

void
DataContainer::
generateInflow()
{
    auto flowValues = parseInflow();

    linearInterpolation(flowValues, M_inflow);
}

void
DataContainer::
linearInterpolation(const std::vector<std::pair<double,double>>& values,
                    std::function<double(double)>& funct)
{
    funct = [values](double x)
    {
        unsigned int count = 0;
        for (auto curpair : values)
        {
            if (x - curpair.first < 1e-15)
                break;
            count++;
        }
        if (count == values.size())
        {
            printlog(YELLOW, "Warning: exiting the bounds of the inflow file");
            return values[count].second;
        }

        if (std::abs(values[count].first - x) < 1e-15)
            return values[count].second;

        double coeff = (values[count].second - values[count-1].second) /
                       (values[count].first - values[count-1].first);

        return values[count-1].second + coeff * (x - values[count-1].first);
    };
}

std::vector<std::pair<double, double>>
DataContainer::
parseInflow()
{
    std::ifstream inflowfile((*M_datafile)("bc_conditions/inflowfile",
                                           "datafiles/inflow.txt"));
    std::vector<std::pair<double, double>> flowValues;
    if (inflowfile.is_open())
    {
        while (inflowfile.good())
        {
            std::pair<double, double> newpair;
            std::string line;
            std::getline(inflowfile,line);
            if (std::strcmp(line.c_str(),""))
            {
                auto firstspace = line.find(" ");
                if (firstspace == std::string::npos || firstspace == 0)
                    throw new Exception("Inflow file is badly formatted");
                newpair.first = std::stod(line.substr(0,firstspace));

                auto lastspace = line.find_last_of(" ");
                if (lastspace == std::string::npos || firstspace == 0)
                    throw new Exception("Inflow file is badly formatted");

                newpair.second = std::stod(line.substr(lastspace+1,line.size()));
                flowValues.push_back(newpair);
            }
        }
        inflowfile.close();
    }
    return flowValues;
}

}
