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

#ifndef BDF_HPP
#define BDF_HPP

#include <redma/RedMA.hpp>
#include <redma/solver/time_marching_algorithms/aTimeMarchingAlgorithm.hpp>
#include <redma/utils/Exception.hpp>

#include <memory>

namespace RedMA
{

template <class DataType>
class BDF : public aTimeMarchingAlgorithm<DataType>
{
public:
    BDF(const GetPot& datafile);

    virtual DataType advance(const double& time, double& dt,
                             SHP(aAssembler) assembler);

    void setup();

protected:
    std::vector<DataType>                   M_prevSolutions;
    std::vector<double>                     M_coefficients;
    unsigned int                            M_order;
    double                                  M_rhsCoeff;
};

}

#include "BDF_imp.hpp"

#endif // BDF_HPP
