//
// Created by Federico Betti on 08/04/2022.
//


#include "LatinHypercube.hpp"

namespace RedMa
{

LatinHypercube::
LatinHypercube(unsigned int numSamples)
{
    M_numSamples = numSamples;

    const double maxAngle = 0.4;
    const double maxAmplitude = 0.2;
    const double maxWidth = 0.3;
    const double minFlow = 0.2;
    const double maxFlow = 0.8;

    M_parametersHandler.registerParameter("in1_alphax", 0,
                                          -maxAngle, maxAngle, true, false);
    M_parametersHandler.registerParameter("in1_alphay", 0,
                                          -maxAngle, maxAngle, true, false);
    M_parametersHandler.registerParameter("in1_alphaz", 0,
                                          -maxAngle, maxAngle, true, false);
    M_parametersHandler.registerParameter("in2_alphax", 0,
                                          -maxAngle, maxAngle, true, false);
    M_parametersHandler.registerParameter("in2_alphay", 0,
                                          -maxAngle, maxAngle, true, false);
    M_parametersHandler.registerParameter("in2_alphaz", 0,
                                          -maxAngle, maxAngle, true, false);
    M_parametersHandler.registerParameter("stenosis_amplitude", 0,
                                          0, maxAmplitude, true, false);
    M_parametersHandler.registerParameter("stenosis_width", 0,
                                          0, maxWidth, true, false);
    M_parametersHandler.registerParameter("flow_rate", 0,
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
generateSamples(unsigned int N, unsigned int d) {
    std::vector<std::map<std::string, double>> LHS_Samples;
    std::vector<std::vector<double>> uniformSamples = drawUniformSamples(N, d);
    std::vector<std::vector<unsigned int>> permutations = getPermutations(N, d);
    for (unsigned int i = 0; i < N; ++i)
    {
        std::vector<double> currentLHS_Sample;
        std::map<std::string, double> currentSample;
        for (unsigned int j = 0; j < d; ++j)
        {
            currentLHS_Sample[j] = 1 / double(N) * (permutations[j][i] + uniformSamples[i][j]);
            currentLHS_Sample[j] = M_paramsBounds[j][0] + currentLHS_Sample[j] * (M_paramsBounds[j][1] - M_paramsBounds[j][0]);
            currentSample.insert(std::pair<std::string, double>(M_paramsNames[j], currentLHS_Sample[j]));
        }
        LHS_Samples.push_back(currentSample);
    }
    std::sort(LHS_Samples.begin(), LHS_Samples.end(), [](const std::map<std::string, double>& sample1,
            const std::map<std::string, double>& sample2)->bool
            { return sample1.at("stenosis_amplitude") < sample2.at("stenosis_amplitude"); });
    return LHS_Samples;
}

std::vector<std::vector<double>>
LatinHypercube::
drawUniformSamples(unsigned int N, unsigned int d) {
    std::vector<std::vector<double>> uniformSamples;
    for (unsigned int i = 0; i < N; ++i)
    {
        std::vector<double> currentSample;
        for (unsigned int j = 0; j < d; ++j)
        {
            std::random_device                  rand_dev;
            std::mt19937                        generator(rand_dev());
            std::uniform_real_distribution<double>  distr(0, 1);
            double randomNumber = distr(generator);
            currentSample.push_back(randomNumber);
        }
        uniformSamples.push_back(currentSample);
    }
    return uniformSamples;
}

std::vector<std::vector<unsigned int>>
LatinHypercube::
getPermutations(unsigned int N, unsigned int d)
{
    std::vector<std::vector<unsigned int>> permutationMatrix;
    std::vector<unsigned int> toBePermuted;
    for (unsigned int i = 0; i < N; ++i) {
        toBePermuted.push_back(i);
    }
    for (unsigned int j = 0; j < d; ++j) {
        std::random_shuffle(toBePermuted.begin(), toBePermuted.end());
        permutationMatrix.push_back(toBePermuted);
    }
    return permutationMatrix;
}
}