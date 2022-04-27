//
// Created by Federico Betti on 08/04/2022.
//

#include "LatinHypercube.hpp"

namespace RedMA
{

LatinHypercube::
LatinHypercube(unsigned int numSamples)
{
    M_numSamples = numSamples;

    const double maxAngle = 0.4;
    const double maxAmplitude = 0.3;
    const double maxWidth = 0.2;
    const double minFlow = 0.2;
    const double maxFlow = 0.8;

    GeometricParametersHandler& parametersHandler = getGeometricParametersHandler();

    parametersHandler.registerParameter("in1_alphax", 0,
                                        -maxAngle, maxAngle, true, false);
    parametersHandler.registerParameter("in1_alphay", 0,
                                        -maxAngle, maxAngle, true, false);
    parametersHandler.registerParameter("in1_alphaz", 0,
                                        -maxAngle, maxAngle, true, false);
    parametersHandler.registerParameter("in2_alphax", 0,
                                        -maxAngle, maxAngle, true, false);
    parametersHandler.registerParameter("in2_alphay", 0,
                                        -maxAngle, maxAngle, true, false);
    parametersHandler.registerParameter("in2_alphaz", 0,
                                        -maxAngle, maxAngle, true, false);
    parametersHandler.registerParameter("stenosis_amplitude", 0,
                                        0, maxAmplitude, true, false);
    parametersHandler.registerParameter("stenosis_width", 0,
                                        0, maxWidth, true, false);
    parametersHandler.registerParameter("flow_rate", 0,
                                        minFlow, maxFlow, true, false);

    setParametersToBeSampled();
}

void
LatinHypercube::
setParametersToBeSampled()
{
    for (auto it = M_parametersHandler.getParametersMap().begin();
    it != M_parametersHandler.getParametersMap().end(); ++ it) {
        M_paramsNames.push_back(it->first);
        M_paramsBounds.push_back(getBounds(it->first));
    }
}

std::vector<double>
LatinHypercube::
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

std::vector<std::map<std::string, double>>
LatinHypercube::
generateSamples(const unsigned int& N, const unsigned int& d)
{
    std::vector<std::map<std::string, double>> LHS_Samples;
    std::vector<std::vector<unsigned int>> permutations = getPermutations(N, d);
    for (unsigned int i = 0; i < N; ++i)
    {
        std::map<std::string, double> currentSample;
        std::vector<double> uniformSample = drawUniformSample(d);
        for (unsigned int j = 0; j < d; ++j)
        {
            double elem = 1 / double(N) * (permutations[j][i] + uniformSample[j]);
            elem = M_paramsBounds[j][0] + elem * (M_paramsBounds[j][1] - M_paramsBounds[j][0]);
            currentSample.insert(std::pair<std::string, double>(M_paramsNames[j], elem));
        }
        LHS_Samples.push_back(currentSample);
    }
    std::sort(LHS_Samples.begin(), LHS_Samples.end(), [](const std::map<std::string, double>& sample1,
            const std::map<std::string, double>& sample2)->bool
            { return sample1.at("stenosis_amplitude") < sample2.at("stenosis_amplitude"); });
    return LHS_Samples;
}

std::vector<double>
LatinHypercube::
drawUniformSample(const unsigned int& d)
{
    std::vector<double> currentSample(d);
    generate(currentSample.begin(), currentSample.end(), []() {return (double) rand()/RAND_MAX; });
    return currentSample;
}

std::vector<std::vector<unsigned int>>
LatinHypercube::
getPermutations(const unsigned int& N, const unsigned int& d)
{
    std::vector<std::vector<unsigned int>> permutationMatrix;
    std::vector<unsigned int> toBePermuted(N);
    std::iota(toBePermuted.begin(), toBePermuted.end(), 0);
    for (unsigned int j = 0; j < d; ++j)
    {
        std::random_shuffle(toBePermuted.begin(), toBePermuted.end());
        permutationMatrix.push_back(toBePermuted);
    }
    return permutationMatrix;
}
}