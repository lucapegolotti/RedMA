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

#ifndef DATACONTAINER_HPP
#define DATACONTAINER_HPP

#include <redma/RedMA.hpp>
#include <redma/utils/Exception.hpp>

#include <lifev/core/filter/GetPot.hpp>

#include <fstream>
#include <functional>
#include <map>

namespace RedMA
{

/*! \brief Container for the problem data.
 *
 * This class is a wrapper for the GetPot containing the problem settings. Moreover,
 * DataContainer takes care of the boundary conditions at the inlet and outlets.
 */
class DataContainer
{
    typedef std::function<double(double)>           Law;
public:
    /// Default constructor
    DataContainer();

    /*! \brief Setter for the datafile.
     *
     * Internally, this method opens the GetPot corresponding to the string
     * passed as argument
     *
     * \param datafile Name of the datafile to be opened (usually, "data").
     */
    void setDatafile(const std::string& datafile);

    /*! \brief  Setter for the inflow law with a given inlet flag
     *
     * \param inflow The inflow rate law.
     * \param inletIndex. The index of the inlet to be considered. If not specified, it defaults to 99 (unique inlet)
     */
    void setInflow(const Law& inflow,
                   unsigned int indexInlet = 99);

    /*! \brief Setter for the distal pressure.
     *
     * \param pressure Analytical function of the distal pressure.
     * \param indexOutlet Index of the outlet to which the distal pressure must
     *        be applied.
     */
    void setDistalPressure(const Law& pressure,
                           const unsigned int& indexOutlet);

    /// Setter for the verbose variable.
    void setVerbose(bool verbose);

    /*! \brief Finalize method.
     *
     * This method must be called whenever the inflow is passed through file.
     * It calls the methods generateRamp() and generateInflow().
     */
    void finalize();

    /// Getter for the internal datafile.
    inline GetPot getDatafile() const {return *M_datafile;}

    /*! \brief Getter for the inflow function corresponding to a specific flag.
     *
     * The inflow function are either set in analytical form through setInflow
     * or passed through file specified in M_datafile.
     *
     * \param flag The flag of the inlet.
     * \return Standard map with key = flag, value = inflow function.
     */
    inline Law getInflow(unsigned int flag = 0) {return M_inflows[flag];}

    /*! \brief Getter for the inflow functions.
     *
     * The inflow function is either set in analytical form through setInflow.
     * or passed through file specified in M_datafile.
     */
    inline std::map<unsigned int, Law> getInflows() {return M_inflows;}

    /// Getter for M_verbose.
    inline bool getVerbose() const {return M_verbose;}

    /*! Getter for the distal pressure at a specific outlet.
     *
     * \param outletIndex Index of the desired outlet.
     */
    std::function<double(double)> getDistalPressure(const unsigned int& outletIndex) const;

    /*! Method to retrieve a char* from the internal datafile.
     *
     * \param location Location within the datafile.
     * \param defValue Default value, returned if location is not found.
     */
    std::string operator()(std::string location, const char* defValue) const;

    /*! Method to retrieve a string from the internal datafile.
     *
     * \param location Location within the datafile.
     * \param defValue Default value, returned if location is not found.
     */
    std::string operator()(std::string location, std::string defValue) const;

    /*! Method to retrieve an int from the internal datafile.
     *
     * \param location Location within the datafile.
     * \param defValue Default value, returned if location is not found.
     */
    int operator()(std::string location, int defValue) const;

    /*! Method to retrieve a double from the internal datafile.
     *
     * \param location Location within the datafile.
     * \param defValue Default value, returned if location is not found.
     */
    double operator()(std::string location, double defValue) const;

    /*! Method to retrieve a bool from the internal datafile.
     *
     * \param location Location within the datafile.
     * \param defValue Default value, returned if location is not found.
     */
    bool operator()(std::string location, bool defValue) const;

    // we differentiate the methods that follow to avoid implicit conversion by +
    // mistake

    /*! Method to set a string in the internal datafile.
     *
     * \param location Location within the datafile.
     * \param value Value to set.
     */
    void setValueString(std::string location, std::string value);

    /*! Method to set an int in the internal datafile.
     *
     * \param location Location within the datafile.
     * \param value Value to set.
     */
    void setValueInt(std::string location, int value);

    /*! Method to set a double in the internal datafile.
     *
     * \param location Location within the datafile.
     * \param value Value to set.
     */
    void setValueDouble(std::string location, double value);

    /*! Method to set a bool in the internal datafile.
     *
     * \param location Location within the datafile.
     * \param value Value to set.
     */
    void setValueBool(std::string location, bool value);

    /*! \brief Evaluate ramp function at specific time.
     *
     * \param time Time in which the ramp must be evaluated.
     */
    double evaluateRamp(double time);

protected:
    std::vector<std::pair<double,double>> parseInflow(std::string filename);

    void generateInflow(std::string inputfilename, unsigned int indexInlet = 99);

    /*! \brief Generate ramp based on the values specified in the time_discretization
     * section of M_datafile.
     */
    void generateRamp();

    /*! \brief Linear interpolation method (used within generateInflow).
     *
     * Given a vector of pairs, this method generate a functional to evaluate
     * the linear interpolant across those coordinates at any location.
     *
     * \param values Pairs of (x,y) values.
     * \funct funct Output param: linear interpolant.
     */
    void linearInterpolation(const std::vector<std::pair<double,double>>& values,
                             std::function<double(double)>& funct);

    /*! \brief Check if the inflow for given inlet should be generated from file or externally set
     *
     * Reads 'generate_inflow' field from data file; if 0 or 1, it returns; if -1
     * (default value), it returns 1 if a datafile for the inflow is available,
     * 0 otherwise.
     *
     * \param indexInlet Index of the inlet to be considered. If not specified, it defaults to 99.
     * \return True if the inflow can be generated from file, false otherwise
     *
     */
    bool checkGenerateInflow(unsigned int indexInlet=99) const;

    shp<GetPot>                                           M_datafile;
    std::map<unsigned int, Law>                           M_inflows;
    Law                                                   M_ramp;
    std::map<unsigned int, Law>                           M_distalPressures;
    bool                                                  M_verbose;
};

}

#endif // DATACONTAINER_HPP
