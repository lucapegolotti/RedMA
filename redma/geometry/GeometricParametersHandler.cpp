#include "GeometricParametersHandler.hpp"

namespace RedMA
{

GeometricParameter::
GeometricParameter(std::string name, const double& value,
                   const double& minValue, const double& maxValue,
                   bool randomizible, bool periodic) :
  M_name(name),
  M_value(value),
  M_originalValue(value),
  M_minValue(minValue),
  M_maxValue(maxValue),
  M_randomizible(randomizible),
  M_periodic(periodic)
{
}

GeometricParameter::
GeometricParameter(const GeometricParameter& other)
{
    M_name = other.M_name;
    M_value = other.M_value;
    M_originalValue = other.M_originalValue;
    M_minValue = other.M_minValue;
    M_maxValue = other.M_maxValue;
    M_randomizible = other.M_randomizible;
    M_periodic = other.M_periodic;
}

std::string
GeometricParameter::
name()
{
    return M_name;
}

int
GeometricParameter::
operator=(const double& value)
{
    if (!M_periodic)
    {
        M_value = value < M_minValue ? M_minValue : value;
        M_value = M_value > M_maxValue ? M_maxValue : M_value;
        if (M_value != value)
        {
            std::string msg = "Parameter with name ";
            msg += M_name;
            msg += " was assigned a value outside the bounds; value = ";
            msg += std::to_string(value);
            msg += ", assigned value = ";
            msg += std::to_string(M_value);
            msg += "\n";
            printlog(RED, msg, true);
            return 1;
        }
    }
    else
    {
        M_value = value;
        double interval = M_maxValue - M_minValue;
        while (M_value > M_maxValue)
            M_value -= interval;

        while (M_value < M_minValue)
            M_value += interval;
    }
    M_originalValue = M_value;
    return 0;
}

void
GeometricParameter::
randomSampleAroundOriginalValue(const double& bounds)
{
    using namespace std;
    // static std::default_random_engine e;
    // static std::uniform_real_distribution<> dis(0, 1);
    // static std::default_random_engine e;
    // e.seed(std::chrono::system_clock::now().time_since_epoch().count());
    // static std::uniform_real_distribution<> dis(0, 1);
    // float randomNumber = dis(e);

    // float randomNumber = dis(e);
    std::random_device rd;
    std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dis(0, 1);
    float randomNumber = dis(gen);
    double minvalue = M_value - bounds;
    double maxvalue = M_value + bounds;

    minvalue = minvalue > M_minValue ? minvalue : M_minValue;
    maxvalue = maxvalue < M_maxValue ? maxvalue : M_maxValue;

    M_value = minvalue + randomNumber * (maxvalue - minvalue);
}

double
GeometricParameter::
getValue()
{
    return M_value;
}

bool
GeometricParameter::
isRandomizible()
{
    return M_randomizible;
}

void
GeometricParameter::
randomSample()
{
    static std::default_random_engine e;
    static std::uniform_real_distribution<> dis(0, 1);

    float randomNumber = dis(e);

    M_value = M_minValue + randomNumber * (M_maxValue - M_minValue);
}

constexpr double GeometricParametersHandler::infty;

GeometricParametersHandler::
GeometricParametersHandler()
{

}

void
GeometricParametersHandler::
registerParameter(std::string name, const double& value,
                  const double& minValue, const double& maxValue,
                  bool randomizible, bool periodic)
{
    if (!exists(name))
    {
        GeometricParameterPtr
          newParameter(new GeometricParameter(name, value, minValue, maxValue,
                                              randomizible, periodic));
        M_parametersMap[name] = newParameter;
    }
    else
    {
        std::string errorMsg = std::string("Cannot register parameter with name ")
                               + name + " multiple times!";
        throw Exception(errorMsg);
    }
}

int
GeometricParametersHandler::
setParameterValue(std::string name, const double& value)
{
    if (exists(name))
    {
        return *M_parametersMap[name] = value;
    }
    else
    {
        std::string errorMsg = std::string("Parameter with name ")
                               + name + " does not exist!";
        throw Exception(errorMsg);
    }
    return 1;
}

void
GeometricParametersHandler::
setRandomValueParameter(std::string name)
{
    M_parametersMap[name]->randomSample();
}

double
GeometricParametersHandler::
operator[](std::string name)
{
    if (!exists(name))
    {
        std::string errorMsg = std::string("Parameter with name ")
                               + name + " does not exist!";
        throw Exception(errorMsg);
    }

    return M_parametersMap[name]->getValue();
}

bool
GeometricParametersHandler::
exists(std::string name)
{
    return M_parametersMap.find(name) != M_parametersMap.end();
}

std::map<std::string,GeometricParametersHandler::GeometricParameterPtr>&
GeometricParametersHandler::
getParametersMap()
{
    return M_parametersMap;
}

void
GeometricParametersHandler::
randomizeParameters()
{
    typedef std::map<std::string, GeometricParameterPtr> mapType;

    for (mapType::iterator it = M_parametersMap.begin();
         it != M_parametersMap.end(); it++)
    {
        GeometricParameterPtr gp = it->second;
        if (gp->isRandomizible())
            gp->randomSample();
    }
}

void
GeometricParametersHandler::
randomizeParametersAroundOriginalValue(const double& bounds)
{
    typedef std::map<std::string, GeometricParameterPtr> mapType;

    for (mapType::iterator it = M_parametersMap.begin();
         it != M_parametersMap.end(); it++)
    {
        GeometricParameterPtr gp = it->second;
        if (gp->isRandomizible())
            gp->randomSampleAroundOriginalValue(bounds);
    }
}

std::vector<double>
GeometricParametersHandler::
getRandomizibleParametersValueAsVector()
{
    typedef std::map<std::string, GeometricParameterPtr> mapType;

    std::vector<double> retVec;

    for (mapType::iterator it = M_parametersMap.begin();
         it != M_parametersMap.end(); it++)
    {
        GeometricParameterPtr gp = it->second;
        if (gp->isRandomizible())
            retVec.push_back(gp->getValue());
    }

    return retVec;
}

}  // namespace RedMA
