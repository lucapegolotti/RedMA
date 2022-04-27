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
#include <chrono>
#include <iostream>
#include <redma/utils/Exception.hpp>
#include <functional>
#include <random>
#include <chrono>

#include <redma/utils/PrintLog.hpp>

#include <ctime>
#include <cstdlib>
#include <redma/RedMA.hpp>

namespace RedMA
{

/// Class for handling a geometric parameter.
class GeometricParameter
{
public:
    /*! Constructor.
     *
     * \param name The name of the geometric parameter.
     * \param value The value of the geometric parameter.
     * \param minValue Lower bound.
     * \param maxValue Upper bound.
     * \param randomizible Set if this parameter can be randomized.
     * \param period Set the parameter to periodic.
     */
    GeometricParameter(std::string name,
                       const double& value,
                       const double& minValue,
                       const double& maxValue,
                       bool randomizible = false,
                       bool periodic = false);

    /*! \brief Copy constructor.
     *
     * \param other The other geometric parameter.
     */
    GeometricParameter(const GeometricParameter& other);

    /*! \brief Equal operator.
     *
     * \param value The value of the parameter.
     * \return Status code; 0 if successful.
     */
    int operator=(const double& value);

    /*! \brief Get value of the parameter.
     *
     * \return The value of the parameter.
     */
    double getValue();

    /// Sample the parameter within the interval specified in the constructor.
    void randomSample();

    /*! \brief Random sample around the value set during construction.
     *
     * \param bounds Bounds around the original value.
     */
    void randomSampleAroundOriginalValue(const double& bounds);

    /*! \brief Return if the parameter is randomizible.
     *
     * \return True if the parameter is randomizible.
     */
    bool isRandomizible();

    /*! \brief Getter for the name of the geometrical parameter.
     *
     * \return The name of the parameter.
     */
    std::string name();

    /*! \brief Getter for the lower bound of this parameter.
     *
     * \return The lower bound of this parameter.
     */
    double getMinValue(){return M_minValue;};

    /*! \brief Getter for the upper bound of this parameter.
     *
     * \return The upper bound of this parameter.
     */
    double getMaxValue(){return M_maxValue;};


private:
    GeometricParameter() {}

    std::string         M_name;
    double              M_value;
    // this is for when we sample around the original value. It is always equal
    // to value, except when we invoke randomSampleAroundOriginalValue
    double              M_originalValue;
    double              M_minValue;
    double              M_maxValue;
    bool                M_randomizible;
    // if the bounds are periodic (e.g. for angles)
    bool                M_periodic;
};

/// Class for handling multiple geometric parameters.
class GeometricParametersHandler
{
    typedef shp<GeometricParameter> GeometricParameterPtr;
public:
    /// Default constructor.
    GeometricParametersHandler();

    /*! \brief Add a new parameter.
     *
     * \param name Name of the new parameter.
     * \param value The value of the parameter.
     * \param minValue The lower bound of the parameter.
     * \param maxValue The upper bound of the parameter.
     * \param randomizible Set the parameter to randomizible.
     * \param periodic Set the parameter to periodic.
     */
    void registerParameter(std::string name,
                           const double& value,
                           const double& minValue,
                           const double& maxValue,
                           bool randomizible = false,
                           bool periodic = false);

    /*! \brief Set the value of a parameter.
     *
     * \param name The name of the parameter.
     * \param value The value of the parameter.
     */
    int setParameterValue(std::string name,
                          const double& value);

    /*! \brief Set a parameter to a random value.
     *
     * \param name The name of the parameter.
     */
    void setRandomValueParameter(std::string name);

    /*! \brief Access operator.
     *
     * \param The name of the parameter we want to access.
     * \return The value of the parameter.
     */
    double operator[](std::string name);

    /*! \brief True if a given parameter exists.
     *
     * \return True if a given parameter exists.
     */
    bool exists(std::string name);

    /*! \brief Get the map of all parameters.
     *
     * \return Parameters map (key = name, value = pointer to the GeometricParameter)
     */
    std::map<std::string,GeometricParameterPtr>& getParametersMap();

    /// Set all randomizible parameters to random values.
    void randomizeParameters();

    /*! \brief Set all randomizible parameters to random values around the original ones.
     *
     *  \param bounds The bound around the original values.
     */
    void randomizeParametersAroundOriginalValue(const double& bounds);

    /*! \brief Getter for a specific parameter.
     *
     * \param name The name of the parameter.
     * \return Shared pointer to the GeometricParameter.
     */
    GeometricParameterPtr getParameter(std::string name){return M_parametersMap[name];};

    /*! \brief Getter for a vector of values of all randomizible parameters.
     *
     * \return The desired vector.
     */
    std::vector<double> getRandomizibleParametersValueAsVector();

    void setParametersFromSample(const std::map<std::string, double>& sample);

    static constexpr double infty = 1e9;

private:
    std::map<std::string, GeometricParameterPtr>    M_parametersMap;
};

}  // namespace RedMA

#endif  // GEOMETRICPARAMETERSHANDLER_HPP
