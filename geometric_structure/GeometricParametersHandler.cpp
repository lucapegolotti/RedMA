#include <GeometricParametersHandler.hpp>

namespace RedMA
{

GeometricParameter::
GeometricParameter(std::string name, const double& value,
                   const double& minValue, const double& maxValue) :
  M_name(name),
  M_value(value),
  M_minValue(minValue),
  M_maxValue(maxValue)
{
}

GeometricParameter::
GeometricParameter(const GeometricParameter& other)
{
    M_name = other.M_name;
    M_value = other.M_value;
    M_minValue = other.M_minValue;
    M_maxValue = other.M_maxValue;
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

void
GeometricParameter::
randomSample()
{
    srand(time(NULL));
    M_value = M_minValue + rand() * (M_maxValue - M_minValue);
}

constexpr double GeometricParametersHandler::infty;

GeometricParametersHandler::
GeometricParametersHandler()
{

}

void
GeometricParametersHandler::
registerParameter(std::string name, const double& value,
                  const double& minValue, const double& maxValue)
{
    if (!exists(name))
    {
        GeometricParameterPtr
          newParameter(new GeometricParameter(name, value, minValue, maxValue));
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

}  // namespace RedMA
