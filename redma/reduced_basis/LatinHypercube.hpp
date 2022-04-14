//
// Created by Federico Betti on 08/04/2022.
//

#ifndef LATINHYPERCUBE_HPP
#define LATINHYPERCUBE_HPP

#include <redma/geometry/GeometricParametersHandler.hpp>
#include <redma/RedMA.hpp>

#include <stdio.h>
#include <cmath>

namespace RedMA
{
    class LatinHypercube
    {
    public:
        LatinHypercube(unsigned int numSamples);

        void setParametersToBeSampled();

        std::vector<double> getBounds(std::string paramName);

        std::vector<std::map<std::string, double>> generateSamples(const unsigned int& N, const unsigned int& d);

        std::vector<double> drawUniformSample(const unsigned int& d);

        std::vector<std::vector<unsigned int>> getPermutations(const unsigned int& N, const unsigned int& d);

        unsigned int getNumSamples() { return M_numSamples; };

        std::vector<std::string> getParamsNames() { return M_paramsNames; };

        unsigned int getNumParams() { return M_paramsNames.size(); };

        std::vector<std::vector<double>> getParamsBounds() { return M_paramsBounds; };

        GeometricParametersHandler& getGeometricParametersHandler() {return M_parametersHandler; };

    private:
        GeometricParametersHandler M_parametersHandler;
        unsigned int M_numSamples;
        std::vector<std::string> M_paramsNames;
        std::vector<std::vector<double>> M_paramsBounds;
    };
} // namespace RedMa

#endif //LATINHYPERCUBE_HPP
