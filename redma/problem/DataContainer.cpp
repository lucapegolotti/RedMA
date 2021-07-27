#include "DataContainer.hpp"

namespace RedMA
{

DataContainer::
DataContainer():
M_verbose(false) {}

void
DataContainer::
setDatafile(const std::string& datafile)
{
    M_datafile.reset(new GetPot(datafile));
}

void
DataContainer::
setInflow(const std::function<double(double)>& inflow,
          unsigned int flag)
{
    bool generate_inflow = this->checkGenerateInflow();

    if (!(generate_inflow))
        M_inflows[flag] = inflow;
    else
        printlog(YELLOW, "[DataContainer] WARNING: Inflow function will be "
                      "read from file, as the 'generate_inflow' flag in datafile is set to 1\n");
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
setValueString(std::string location, std::string value)
{
    M_datafile->set(location.c_str(), value.c_str());
}

void
DataContainer::
setValueInt(std::string location, int value)
{
    M_datafile->set(location.c_str(), value);
}

void
DataContainer::
setValueDouble(std::string location, double value)
{
    M_datafile->set(location.c_str(), value);
}

void
DataContainer::
setValueBool(std::string location, bool value)
{
    M_datafile->set(location.c_str(), value);
}

void
DataContainer::
finalize()
{
    if (M_inflows.size() == 0)
    {
        generateRamp();

        int ninlets = (*M_datafile)("bc_conditions/numinletbcs", -1);
        if (ninlets == -1)
        {
            std::string inputfile = (*M_datafile)("bc_conditions/inflowfile",
                                         "datafiles/inflow.txt");
            generateInflow(0, inputfile);
        }
        else
        {
            for (unsigned int i = 0; i < ninlets; i++)
            {
                std::string path = "bc_conditions/inlet" + std::to_string(i);
                std::string arg1 = path + "/flag";
                unsigned int curflag = (*M_datafile)(arg1.c_str(), 0);
                arg1 = path + "/inflowfile";
                std::string arg2 = "datafiles/inflow" + std::to_string(i) + ".txt";
                std::string inputfile = (*M_datafile)(arg1.c_str(),arg2.c_str());
                generateInflow(curflag, inputfile);
            }
        }
    }

    if (!M_inflow)
        throw new Exception("An inflow function has neither being set nor being "
                            "interpolated from datafile! Either call to  'setInflow' method before "
                            "the 'finalize' method (with 'generate_inflow' flag set to 0) or call "
                            "the 'finalize' method (with 'generate_inflow' flag set to 1 and "
                            "providing a valid inflow text file)!");
}

void
DataContainer::
generateInflow(unsigned int flag, std::string inputfilename)
{
    bool generate_inflow = this->checkGenerateInflow();

    if (generate_inflow)
    {
        auto flowValues = parseInflow(inputfilename);
        linearInterpolation(flowValues, M_inflows[flag])
    }
}

void
DataContainer::
generateRamp()
{
    double t0ramp = (*M_datafile)("time_discretization/t0ramp", 0.0);
    double t0 = (*M_datafile)("time_discretization/t0", 0.0);

    if (std::abs(t0ramp - t0) > 1e-15)
    {
        M_ramp = [t0,t0ramp](double x)
        {
            return (1.0 - std::cos((x - t0ramp) * M_PI / (t0 - t0ramp))) / 2.0;
        };
    }
}

void
DataContainer::
linearInterpolation(const std::vector<std::pair<double,double>>& values,
                    std::function<double(double)>& funct)
{
    funct = [values,this](double x)
    {
        if (x < values[0].first && M_ramp)
            return M_ramp(x) * values[0].second;

        unsigned int count = 0;
        for (auto curpair : values)
        {
            if (x - curpair.first < 1e-15)
                break;
            count++;
        }
        if (count == values.size())
        {
            printlog(YELLOW, "WARNING: exiting the bounds of the inflow file");
            return values[count-1].second;
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
parseInflow(std::string filename)
{
    std::ifstream inflowfile;
    inflowfile.open(filename);


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

double
DataContainer::
evaluateRamp(double time)
{
    double t0 = (*M_datafile)("time_discretization/t0", 0.0);

    if (time < t0 && M_ramp)
        return M_ramp(time);
    else
        return 1;
}

bool
DataContainer::
checkGenerateInflow() const
{
    int generate_inflow = (*M_datafile)("bc_conditions/generate_inflow", -1);

    if (generate_inflow == -1)
    {
        std::ifstream inflowfile((*M_datafile)("bc_conditions/inflowfile",
                                               "datafiles/inflow.txt"));
        generate_inflow = (inflowfile.good()) ? 1 : 0;
    }

    return (generate_inflow == 1);
}
}
