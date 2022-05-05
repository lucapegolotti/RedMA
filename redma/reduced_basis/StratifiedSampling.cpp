//
// Created by betti on 5/4/22.
//

#include "StratifiedSampling.hpp"

namespace RedMA
{

StratifiedSampling::
StratifiedSampling(unsigned int numSamples)
{
    M_numSamples = numSamples;

    const double maxAngle = 0.4;
    const double maxAmplitude = 0.04;
    const double maxWidth = 0.04;
    const double minFlow = 0.6;
    const double maxFlow = 0.64;

    GeometricParametersHandler& paramsHandler = getParametersHandler();

//    paramsHandler.registerParameter("in1_alphax", 0,
//                                          -maxAngle, maxAngle, true, false);
//    paramsHandler.registerParameter("in1_alphay", 0,
//                                          -maxAngle, maxAngle, true, false);
//    paramsHandler.registerParameter("in1_alphaz", 0,
//                                          -maxAngle, maxAngle, true, false);
//    paramsHandler.registerParameter("in2_alphax", 0,
//                                          -maxAngle, maxAngle, true, false);
//    paramsHandler.registerParameter("in2_alphay", 0,
//                                          -maxAngle, maxAngle, true, false);
//    paramsHandler.registerParameter("in2_alphaz", 0,
//                                          -maxAngle, maxAngle, true, false);

    paramsHandler.registerParameter("flow_rate", 0,
                                    minFlow, maxFlow, true, false);
    paramsHandler.registerParameter("stenosis_width", 0,
                                    0, maxWidth, true, false);
    paramsHandler.registerParameter("stenosis_amplitude", 0,
                                    0, maxAmplitude, true, false);


    setParametersToBeSampled();
}

void
StratifiedSampling::
setParametersToBeSampled()
{
    for (auto it = M_parametersHandler.getParametersMap().begin();
         it != M_parametersHandler.getParametersMap().end(); ++ it) {
        M_paramsNames.push_back(it->first);
        M_paramsBounds.push_back(getBounds(it->first));
    }
}

std::vector<double>
StratifiedSampling::
getBounds(std::string paramName)
{
    std::vector<double> bounds;
    if (M_parametersHandler.exists(paramName))
    {
        bounds.push_back(M_parametersHandler.getParametersMap().at(paramName)->getMinValue());
        bounds.push_back(M_parametersHandler.getParametersMap().at(paramName)->getMaxValue());
    }
    return bounds;
}

std::map<std::string, std::vector<double>>
StratifiedSampling::
generateSamples(unsigned int N, unsigned int d)
{
    std::map<std::string, std::vector<double>> samples;
    for (unsigned int j = 0; j < d; ++j)
    {
        std::vector<double> oneParamSamples;
        for (unsigned int i = 0; i < N; ++i)
        {
            double elem = M_paramsBounds[j][0] + i/float(N) * (M_paramsBounds[j][1] - M_paramsBounds[j][0]);
            oneParamSamples.push_back(elem);
        }
        samples.insert(std::pair<std::string, std::vector<double>> (M_paramsNames[j], oneParamSamples));
    }
    return samples;
}

}
