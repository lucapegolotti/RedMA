//
// Created by Federico Betti on 08/04/2022.
//

#ifndef LATINHYPERCUBE_HPP
#define LATINHYPERCUBE_HPP

#include <redma/RedMA.hpp>

#include <cmath>
#include <iomanip>
#include <fstream>
#include <random>

namespace RedMa {
    class LatinHypercube {
    public:
        LatinHypercube(unsigned int numSamples, unsigned int dim);

        std::vector<std::vector<double>> generateSamples(unsigned int N, unsigned int);

        std::vector<std::vector<double>> drawUniformSamples(unsigned int N, unsigned int d);

        std::vector<std::vector<unsigned int>> getPermutations(unsigned int N, unsigned int d);

        unsigned int getNumSamples() { return M_numSamples; };

        unsigned int getDimension() { return M_dimension; };

    private:
        unsigned int M_numSamples;
        unsigned int M_dimension;
    };
} // namespace RedMa

#endif //LATINHYPERCUBE_HPP
