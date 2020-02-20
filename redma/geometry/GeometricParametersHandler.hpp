// Reduced Modeling of Arteries (RedMA)
// Copyright (C) 2019  Luca Pegolotti
//
// RedMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RedMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef GEOMETRICPARAMETERSHANDLER_HPP
#define GEOMETRICPARAMETERSHANDLER_HPP

#include <map>
#include <memory>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <redma/utils/Exception.hpp>

#include <redma/utils/PrintLog.hpp>

#include <ctime>
#include <cstdlib>

namespace RedMA
{

class GeometricParameter
{
public:
    GeometricParameter(std::string name, const double& value,
                       const double& minValue, const double& maxValue,
                       bool randomizible = false,
                       bool periodic = false);

    GeometricParameter(const GeometricParameter& other);

    int operator=(const double& value);

    double getValue();

    void randomSample();

    void randomSampleAroundOriginalValue(const double& bounds);

    bool isRandomizible();

    std::string name();

    double getMinValue(){return M_minValue;};

    double getMaxValue(){return M_maxValue;};

private:
    GeometricParameter() {}

    std::string M_name;
    double M_value;
    // this is for when we sample around the original value. It is always equal
    // to value, except when we invoke randomSampleAroundOriginalValue
    double M_originalValue;
    double M_minValue;
    double M_maxValue;
    bool M_randomizible;
    // if the bounds are periodic (e.g. for angles)
    bool M_periodic;
};

class GeometricParametersHandler
{
    typedef std::shared_ptr<GeometricParameter> GeometricParameterPtr;
public:
    GeometricParametersHandler();

    void registerParameter(std::string name, const double& value,
                           const double& minValue, const double& maxValue,
                           bool randomizible = false,
                           bool period = false);

    int setParameterValue(std::string name, const double& value);

    void setRandomValueParameter(std::string name);

    double operator[](std::string name);

    bool exists(std::string name);

    static constexpr double infty = 1e9;

    std::map<std::string,GeometricParameterPtr>& getParametersMap();

    void randomizeParameters();

    void randomizeParametersAroundOriginalValue(const double& bounds);

    GeometricParameterPtr getParameter(std::string name){return M_parametersMap[name];};

private:
    std::map<std::string, GeometricParameterPtr> M_parametersMap;
};

}  // namespace RedMA

#endif  // GEOMETRICPARAMETERSHANDLER_HPP
