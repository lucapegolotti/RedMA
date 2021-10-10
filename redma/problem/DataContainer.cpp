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
setInletBC(const std::function<double(double)>& inletLaw,
           unsigned int indexInlet)
{
    bool generate_inletBC = this->checkGenerateInletBC(indexInlet);

    if (!(generate_inletBC)) {
        std::string path;
        if (indexInlet != 99)
            path = "bc_conditions/inlet" + std::to_string(indexInlet) +"/flag";
        else
            path = "bc_conditions/flag";

        unsigned int flag = (*M_datafile)(path.c_str(), 0);
        M_inletBCs[flag] = inletLaw;
    }
    else
        printlog(YELLOW, "[DataContainer] WARNING: Inlet BC function will be "
                      "read from file, as the 'generate_inletBC' flag in datafile is set to 1\n");
}

void
DataContainer::
setOutletBC(const std::function<double(double)>& outletLaw,
           unsigned int indexOutlet)
{
    M_outletBCs[indexOutlet] = outletLaw;
}

DataContainer::Law
DataContainer::
getDistalPressure(const unsigned int& outletIndex) const
{
    auto it = M_distalPressures.find(outletIndex);

    if (it == M_distalPressures.end())
    {
        printlog(YELLOW, "[DataContainer] WARNING: distal pressure not set in outlet number " +
                std::to_string(outletIndex) + ". Setting a default null distal pressure.\n");
        std::function<double(double)> nullP = [](double t) {return 0.0;};
        return nullP;
    }
    else
    {
        return it->second;
    }
}

DataContainer::Law
DataContainer::
getIntramyocardialPressure() const
{
    return M_intraMyocardialPressure;
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
    // handling inlets
    if (M_inletBCs.empty())
    {
        generateRamp();
        std::string BC_type = (*M_datafile)("bc_conditions/inlet_bc_type", "dirichlet");

        int ninlets = (*M_datafile)("bc_conditions/numinletbcs", -1);
        if (ninlets == -1)
        {
            std::string inputfile;
            if (!(std::strcmp(BC_type.c_str(), "dirichlet")))
                inputfile = (*M_datafile)("bc_conditions/inflowfile",
                            "datafiles/inflow.txt");
            else if (!(std::strcmp(BC_type.c_str(), "neumann")))
                inputfile = (*M_datafile)("bc_conditions/inletpressurefile",
                        "datafiles/inlet_pressure.txt");
            else
                throw new Exception("Unrecognized inlet BC type " + BC_type);

            generateInletBC(inputfile, 99);
        }
        else
        {
            for (unsigned int i = 0; i < ninlets; i++)
            {
                std::string path = "bc_conditions/inlet" + std::to_string(i);
                std::string arg1;
                std::string arg2;
                if (!(std::strcmp(BC_type.c_str(), "dirichlet"))) {
                    arg1 = path + "/inflowfile";
                    arg2 = "datafiles/inflow" + std::to_string(i) + ".txt";
                }
                else if (!(std::strcmp(BC_type.c_str(), "neumann"))) {
                    arg1 = path + "/inletpressurefile";
                    arg2 = "datafiles/inlet_pressure" + std::to_string(i) + ".txt";
                }
                else
                    throw new Exception("Unrecognized inlet BC type " + BC_type);

                std::string inputfile = (*M_datafile)(arg1.c_str(),arg2.c_str());
                generateInletBC(inputfile, i);
            }
        }
    }

    if (M_inletBCs.empty())
        throw new Exception("An inflow function has neither being set nor being "
                            "interpolated from datafile! Either call to  'setInletBC' method before "
                            "the 'finalize' method (with 'generate_inletBC' flag set to 0) or call "
                            "the 'finalize' method (with 'generate_inletBC' flag set to 1 and "
                            "providing a valid inflow text file)!");

    // handling outlets
    unsigned int noutletbcs = (*M_datafile)("bc_conditions/numoutletbcs", 0);
    std::string dataEntry;
    std::string BCType;
    unsigned int indexOutlet = 0;
    while (indexOutlet < noutletbcs)
    {
        dataEntry = "bc_conditions/outlet" + std::to_string(indexOutlet) +"/type";
        BCType = (*M_datafile)(dataEntry.c_str(), "windkessel");
        if (!std::strcmp(BCType.c_str(), "coronary"))
            indexOutlet = noutletbcs + 2;
        else
            indexOutlet++;
    }
    if (indexOutlet == noutletbcs + 2)
    {
        std::string IMPpath = (*M_datafile)("bc_conditions/pimfile", "datafiles/IM_pressure.txt");
        this->generateIntraMyocardialPressure(IMPpath);
    }
}

void
DataContainer::
generateInletBC(std::string inputfilename, unsigned int indexInlet)
{
    bool generate_inletBC = this->checkGenerateInletBC(indexInlet);

    if (generate_inletBC)
    {
        std::string path;
        if (indexInlet != 99)
            path = "bc_conditions/inlet" + std::to_string(indexInlet) +"/flag";
        else
            path = "bc_conditions/flag";
        unsigned int flag = (*M_datafile)(path.c_str(), 1);

        auto values = parseTimeValueFile(inputfilename);
        linearInterpolation(values, M_inletBCs[flag]);

        std::string BC_type = (*M_datafile)("bc_conditions/inlet_bc_type", "dirichlet");
        if (!std::strcmp(BC_type.c_str(), "neumann")) {
            std::function<double(double)> tmp = M_inletBCs[flag];
            M_inletBCs[flag] = [tmp](double x) {return tmp(x) * (-1.0);};
        }
    }
}

void
DataContainer::
generateIntraMyocardialPressure(std::string inputfilename)
{
    try
    {
        auto values = parseTimeValueFile(inputfilename);
        linearInterpolation(values, M_intraMyocardialPressure);
    }
    catch (Exception* e)
    {
        printlog(YELLOW, "[DataContainer] Intramyocardial pressure datafile not found; "
                         "setting it to 0 by default.");
        M_intraMyocardialPressure = [](double t) {return 0.0;};
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
        M_ramp = [t0,t0ramp](double t)
        {
            if (t <= t0)
                return (1.0 - std::cos((t - t0ramp) * M_PI / (t0 - t0ramp))) / 2.0;
            else
                return 1.0;
        };
    }
}

void
DataContainer::
linearInterpolation(const std::vector<std::pair<double,double>>& values,
                    std::function<double(double)>& funct)
{
    if (std::abs(values[0].second) <= 1e-5 && M_ramp)
        printlog(YELLOW, "[DataContainer] WARNING: unnecessary addition of ramp, as the initial "
                         "value is already zero. Risk of numerical problems!\n");

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
            printlog(YELLOW, "[DataContainer] WARNING: exiting the bounds of the inflow file");
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
parseTimeValueFile(std::string filename)
{
    std::ifstream inletfile;
    inletfile.open(filename);

    std::vector<std::pair<double, double>> values;
    if (inletfile.is_open())
    {
        while (inletfile.good())
        {
            std::pair<double, double> newpair;
            std::string line;
            std::getline(inletfile>>std::ws,line);

            if (!(line.empty()))
            {
                std::stringstream lineStream(line);
                lineStream >> newpair.first >> newpair.second;

                if (!inletfile || !lineStream)
                    throw new Exception("Failed in reading file " + filename);

                values.push_back(newpair);
            }
        }
        inletfile.close();
    }
    return values;
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
checkGenerateInletBC(unsigned int indexInlet) const
{
    std::string inlet_path;
    std::string default_path;
    std::string BC_type = (*M_datafile)("bc_conditions/inlet_bc_type", "dirichlet");
    if (indexInlet != 99)
    {
        inlet_path = "bc_conditions/inlet" + std::to_string(indexInlet);
        if (!(std::strcmp(BC_type.c_str(), "dirichlet")))
            default_path = "datafiles/inflow" + std::to_string(indexInlet) + ".txt";
        else if (!(std::strcmp(BC_type.c_str(), "neumann")))
            default_path = "datafiles/inlet_pressure" + std::to_string(indexInlet) + ".txt";
        else
            throw new Exception("Unrecognized inlet BC type " + BC_type);
    }
    else
    {
        inlet_path = "bc_conditions";
        if (!(std::strcmp(BC_type.c_str(), "dirichlet")))
            default_path = "datafiles/inflow.txt";
        else if (!(std::strcmp(BC_type.c_str(), "neumann")))
            default_path = "datafiles/inlet_pressure.txt";
        else
            throw new Exception("Unrecognized inlet BC type " + BC_type);
    }

    std::string path1 = inlet_path + "/generate_inletBC";
    int generate_inletBC = (*M_datafile)(path1.c_str(), -1);

    if (generate_inletBC == -1)
    {
        std::string path2;
        if (!(std::strcmp(BC_type.c_str(), "dirichlet")))
            path2 = inlet_path + "/inflowfile";
        else if (!(std::strcmp(BC_type.c_str(), "neumann")))
            path2 = inlet_path + "/inletpressurefile";
        else
            throw new Exception("Unrecognized inlet BC type " + BC_type);

        std::ifstream inflowfile((*M_datafile)(path2.c_str(),
                                                 default_path.c_str()));
        generate_inletBC = (inflowfile.good()) ? 1 : 0;
    }

    return (generate_inletBC == 1);
}
}
