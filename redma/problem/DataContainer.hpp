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

    /*! \brief  Setter for the inlet law (either Dirichlet or Neumann) with a given inlet flag
     *
     * \param inflow The inlet law.
     * \param inletIndex. The index of the inlet to be considered. If not specified, it defaults to 99 (unique inlet)
     */
    void setInletBC(const Law& inflow,
                    unsigned int indexInlet = 99);

    /*! \brief  Setter for the outlet law (Neumann) with a given outlet index
     *
     * \param inflow The outlet law.
     * \param inletIndex. The index of the outlet to be considered
     */
    void setOutletBC(const Law& inflow,
                    unsigned int indexOutlet);

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
     * It calls the methods generateRamp() and generateInletBC().
     */
    void finalize();

    /// Getter for the internal datafile.
    inline GetPot getDatafile() const {return *M_datafile;}

    /*! \brief Getter for the inlet function corresponding to a specific flag.
     *
     * The inlet function is either set in analytical form through setInletBC
     * or passed through file specified in M_datafile.
     *
     * \param flag The flag of the inlet.
     * \return Inlet function corresponding to the given flag
     */
    inline Law getInletBC(unsigned int flag = 0) {return M_inletBCs[flag];}

    /*! \brief Getter for the inflow functions.
     *
     * The inlet function is either set in analytical form through setInletBC.
     * or passed through file specified in M_datafile.
     *
     * \return Standard map with key = flag, value = inlet function.
     */
    inline std::map<unsigned int, Law> getInletBCs() {return M_inletBCs;}

    /*! \brief Getter for the outlet function corresponding to a specific flag.
     *
     * The outlet function is set in analytical form through setOutletBC
     *
     * \param outletIndex The index of the outlet.
     * \return Outlet function corresponding to the given outlet index
     */
    inline Law getOutletBC(unsigned int outletIndex = 0) {return M_outletBCs[outletIndex];}

    /*! \brief Getter for the outlet functions.
     *
     * The outlet function is set in analytical form through setOutletBC.
     *
     * \return Standard map with key = index, value = outlet function.
     */
    inline std::map<unsigned int, Law> getOutletBCs() {return M_outletBCs;}

    /// Getter for M_verbose.
    inline bool getVerbose() const {return M_verbose;}

    /*! \brief Getter for the distal pressure at a specific outlet.
     *
     * \param outletIndex Index of the desired outlet.
     * \return function representing the distal pressure at the prescribed outlet
     */
    Law getDistalPressure(const unsigned int& outletIndex) const;

    /*! \brief Getter for the intramyocardial pressure
     *
     * \return Function describing the intramyocardial pressure at the prescribed outlet
     */
    Law getIntramyocardialPressure() const;

    /*! \brief Method to retrieve a char* from the internal datafile.
     *
     * \param location Location within the datafile.
     * \param defValue Default value, returned if location is not found.
     */
    std::string operator()(std::string location, const char* defValue) const;

    /*! \brief Method to retrieve a string from the internal datafile.
     *
     * \param location Location within the datafile.
     * \param defValue Default value, returned if location is not found.
     */
    std::string operator()(std::string location, std::string defValue) const;

    /*! \brief Method to retrieve an int from the internal datafile.
     *
     * \param location Location within the datafile.
     * \param defValue Default value, returned if location is not found.
     */
    int operator()(std::string location, int defValue) const;

    /*! \brief Method to retrieve a double from the internal datafile.
     *
     * \param location Location within the datafile.
     * \param defValue Default value, returned if location is not found.
     */
    double operator()(std::string location, double defValue) const;

    /*! \brief Method to retrieve a bool from the internal datafile.
     *
     * \param location Location within the datafile.
     * \param defValue Default value, returned if location is not found.
     */
    bool operator()(std::string location, bool defValue) const;

    // we differentiate the methods that follow to avoid implicit conversion by +
    // mistake

    /*! \brief Method to set a string in the internal datafile.
     *
     * \param location Location within the datafile.
     * \param value Value to set.
     */
    void setValueString(std::string location, std::string value);

    /*! \brief Method to set an int in the internal datafile.
     *
     * \param location Location within the datafile.
     * \param value Value to set.
     */
    void setValueInt(std::string location, int value);

    /*! \brief Method to set a double in the internal datafile.
     *
     * \param location Location within the datafile.
     * \param value Value to set.
     */
    void setValueDouble(std::string location, double value);

    /*! \brief Method to set a bool in the internal datafile.
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

    /*! \brief Parse text file storing pairs of values in the form "time  value"
     *
     * \param filename Name of the file storing the data
     * \return Vector of pairs of the form (time,value)
     */
    std::vector<std::pair<double,double>> parseTimeValueFile(std::string filename);

    /*! \brief Generate law at the specified inlet by parsing data file and linear interpolation
     *
     * \param inputfilename Name of the file storing the inlet BC
     * \param indexInlet Index of the inlet to be consider. If the inlet in unique, it defaults to 99.
     */
    void generateInletBC(std::string inputfilename, unsigned int indexInlet = 99);

    /*! \brief Generate by linear interpolation from file and return the IntraMyocardial pressure.
     *
     * \param inputfilename Name of the file storing the data of Intramyocardial pressure
     */
    void generateIntraMyocardialPressure(std::string inputfilename);

    /*! \brief Generate ramp based on the values specified in the time_discretization
     * section of M_datafile.
     */
    void generateRamp();

    /*! \brief Linear interpolation method (used within generateInletBC).
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
     * Reads 'generate_inletBC' field from data file; if 0 or 1, it returns; if -1
     * (default value), it returns 1 if a datafile for the inflow is available,
     * 0 otherwise.
     *
     * \param indexInlet Index of the inlet to be considered. If not specified, it defaults to 99.
     * \return True if the inflow can be generated from file, false otherwise
     *
     */
    bool checkGenerateInletBC(unsigned int indexInlet=99) const;

    shp<GetPot>                                           M_datafile;
    std::map<unsigned int, Law>                           M_inletBCs;
    std::map<unsigned int, Law>                           M_outletBCs;
    Law                                                   M_ramp;
    std::map<unsigned int, Law>                           M_distalPressures;
    Law                                                   M_intraMyocardialPressure;
    bool                                                  M_verbose;
};

}

#endif // DATACONTAINER_HPP
