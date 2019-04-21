#include <GeometricParametersHandler.hpp>

namespace RedMA
{

GeometricParameter::
GeometricParameter(std::string name, const double& value,
                   const double& minValue, const double& maxValue,
                   bool randomizible) :
  M_name(name),
  M_value(value),
  M_minValue(minValue),
  M_maxValue(maxValue),
  M_randomizible(randomizible)
{
}

GeometricParameter::
GeometricParameter(const GeometricParameter& other)
{
    M_name = other.M_name;
    M_value = other.M_value;
    M_minValue = other.M_minValue;
    M_maxValue = other.M_maxValue;
    M_randomizible = other.M_randomizible;
}

std::string
GeometricParameter::
name()
{
    return M_name;
}

void
GeometricParameter::
operator=(const double& value)
{
    M_value = value < M_minValue ? M_minValue : value;
    M_value = M_value > M_maxValue ? M_maxValue : M_value;
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
    // note: the seed must have been set at this point
    float randomNumber =
                        static_cast<float>(rand())/static_cast<float>(RAND_MAX);
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
                  bool randomizible)
{
    if (!exists(name))
    {
        GeometricParameterPtr
          newParameter(new GeometricParameter(name, value, minValue, maxValue,
                                              randomizible));
        M_parametersMap[name] = newParameter;
    }
    else
    {
        std::string errorMsg = std::string("Cannot register parameter with name ")
                               + name + " multiple times!";
        throw Exception(errorMsg);
    }
}

void
GeometricParametersHandler::
setParameterValue(std::string name, const double& value)
{
    if (exists(name))
    {
        *M_parametersMap[name] = value;
    }
    else
    {
        std::string errorMsg = std::string("Parameter with name ")
                               + name + " does not exist!";
        throw Exception(errorMsg);
    }
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

}  // namespace RedMA
