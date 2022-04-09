//
// Created by Federico Betti on 08/04/2022.
//


#include "LatinHypercube.hpp"
#include <redma/RedMA.hpp>

namespace RedMa
{

LatinHypercube::
LatinHypercube(unsigned int numSamples, unsigned int dim) {
    M_numSamples = numSamples;
    M_dimension = dim;
}

std::vector<std::vector<double>>
LatinHypercube::
generateSamples(unsigned int N, unsigned int d) {
    std::vector<std::vector<double>> LHS_Samples;
    std::vector<std::vector<double>> uniformSamples = drawUniformSamples(N, d);
    std::vector<std::vector<unsigned int>> permutations = getPermutations(N, d);
    for (unsigned int i = 0; i < N; ++i) {
        std::vector<double> currentLHS_Sample;
        for (unsigned int j = 0; j < d; ++j)
        {
            currentLHS_Sample[j] = 1 / double(N) * (permutations[j][i] - 1 + uniformSamples[i][j]);
        }
        LHS_Samples.push_back(currentLHS_Sample);
    }

    return LHS_Samples;
}

std::vector<std::vector<double>>
LatinHypercube::
drawUniformSamples(unsigned int N, unsigned int d) {
    std::vector<std::vector<double>> uniformSamples;
    for (unsigned int i = 0; i < N; ++i) {
        std::vector<double> currentSample;
        for (unsigned int j = 0; j < d; ++j) {
            static std::default_random_engine e;
            static std::uniform_real_distribution<> dis(0, 1);

            double randomNumber = dis(e);
            currentSample[j] = randomNumber;
        }
        uniformSamples.push_back(currentSample);
    }
    return uniformSamples;
}

std::vector<std::vector<unsigned int>>
LatinHypercube::
getPermutations(unsigned int N, unsigned int d) {
    std::vector<std::vector<unsigned int>> permutationMatrix;
    std::vector<unsigned int> toBePermuted;
    for (unsigned int i = 0; i < N; ++i) {
        toBePermuted.push_back(i);
    }
    for (unsigned int j = 0; j < d; ++j) {
        std::vector<unsigned int> CopyPermuted = toBePermuted;
        std::random_shuffle(CopyPermuted.begin(), CopyPermuted.end());
        permutationMatrix.push_back(CopyPermuted);
    }
    return permutationMatrix;
}
}