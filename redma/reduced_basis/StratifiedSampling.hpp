//
// Created by betti on 5/4/22.
//

#ifndef STRATIFIEDSAMPLING_HPP
#define STRATIFIEDSAMPLING_HPP

#include <redma/RedMA.hpp>
#include <redma/problem/DataContainer.hpp>
#include <redma/geometry/building_blocks/Bypass.hpp>
#include <redma/assemblers/finite_element/StokesAssemblerFE.hpp>
#include <redma/geometry/GeometricParametersHandler.hpp>

namespace RedMA
{
    class StratifiedSampling {
    public:
        StratifiedSampling(std::vector<unsigned int> numSamples);

        void setParametersToBeSampled();

        std::vector<double> getBounds(std::string paramName);

        std::vector<double> setParametersStepSize();

        std::map<std::string, std::vector<double>> generateSamples();

        unsigned int getNumSamples() { return M_numSamples; };

        std::vector<std::string> getParamsNames() { return M_paramsNames; };

        unsigned int getNumParams() { return M_paramsNames.size(); };

        std::vector<std::vector<double>> getParamsBounds() { return M_paramsBounds; };

        GeometricParametersHandler& getParametersHandler() { return M_parametersHandler; };

        std::vector<unsigned int> getNumComponents() { return M_numPerComponent; };

    private:
        RedMA::GeometricParametersHandler M_parametersHandler;
        unsigned int M_numSamples;
        std::vector<unsigned int> M_numPerComponent;
        std::vector<std::string> M_paramsNames;
        std::vector<std::vector<double>> M_paramsBounds;
    };
}

#endif //STRATIFIEDSAMPLING_HPP
