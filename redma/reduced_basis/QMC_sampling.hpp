//
// Created by Federico Betti on 03/04/2022.
//

#ifndef QMC_SAMPLING_H
#define QMC_SAMPLING_H

#include <redma/geometry/GeometricParametersHandler.hpp>
#include <redma/RedMA.hpp>
#include <redma/geometry/building_blocks/Bypass.hpp>
#include <redma/assemblers/finite_element/StokesAssemblerFE.hpp>
#include <redma/problem/DataContainer.hpp>


// make getters and reset all the attributes as private, then change accordingly in the snapshots class

namespace RedMA
{
    /// sampling class for geometric and physical parameters

    class QMC_sampling
    {
    public:
        /*! \brief Default constructor
        * \param numSamples the number of samples
        */
        QMC_sampling(unsigned int numSamples);

        /*! \brief Returns the bounds for the parameters
         * \param paramName the name of the parameter
         * \param paramsHandler handler of the parameters
         */
        std::vector<double> getBounds(std::string paramName);

        /*! \brief Returns the generating vector of the lattice point set
        */
        std::vector<double> getGeneratingVector() {return M_generatingVector; };

        /*! \brief Returns the parameters bounds
        */
        std::vector<std::vector<double>> getParamsBounds() {return M_paramsBounds; };

        /*! \brief Returns the parameters names
        */
        std::vector<std::string> getParamsNames() {return M_paramsNames; };

        /*! \brief Returns the number of samples
        */
        unsigned int getNumSamples() {return M_numSamples; };

        /*! \brief This function performs the sampling given in input the number of samples desired.
         * \param N the number of samples
        */

        std::vector<std::map<std::string, double>> getSamples(unsigned int N, std::vector<double> generatingVector,
                   std::vector<std::string> paramsNames, std::vector<std::vector<double>> paramsBounds);

        /*! \brief This function performs component-wise multiplication of a vector by a scalar
         *
         * @param v vector
         * @param k scalar
         * @return k * v
         */
        std::vector<double> MultiplyVectorByScalar(std::vector<double> v, double k);

        /*! \brief This function gets component-wise the fractional part of a vector
         *
         * @param v vector
         * @return {v}
         */
        std::vector<double> getFractionalPart(std::vector<double> v);

        /*! \brief This function maps component-wise the samples to the intervals of definition
         *
         * @param v sample
         * @param paramsBounds bounds for each component (parameter)
         * @return a + v*(b-a)
         */
        std::vector<double> translateSampleToBounds(std::vector<double> v, std::vector<std::vector<double>> paramsBounds);

        /*! \brief This function creates the map of the parameters subject to sampling
         *
         * @param paramsNames vector with strings all the names of the parameters
         * @param paramsHandler handler for the bounds
         * @return final vector of the parameters subject to sampling with bounds in a vector of Vector2D
         */
        void setParametersToBeSampled();

        /*! \brief simple checks about consistency of the initializations
         *
         * @param paramsMap bounds map
         * @param generatingVector generating vector
         */
        void checkConsistency();

        /*! \brief getter for the parameters handler
         *
         * @return M_parametersHandler
         */
        GeometricParametersHandler& getGeometricParametersHandler() {return M_parametersHandler; };

    private:
        std::string M_name;
        GeometricParametersHandler M_parametersHandler;
        std::vector<double> M_generatingVector;
        unsigned int M_numSamples;
        std::vector<std::vector<double>> M_paramsBounds;
        std::vector<std::string> M_paramsNames;
        std::vector<std::map<std::string, double>> M_samplesVector;
    };
}

#endif //QMC_SAMPLING_H

