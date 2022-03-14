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

#ifndef APRESSUREDROP_HPP
#define APRESSUREDROP_HPP

#include <redma/RedMA.hpp>

#include <redma/solver/time_marching_algorithms/aFunctionProvider.hpp>
#include <redma/array/DoubleVector.hpp>
#include <redma/array/DoubleMatrix.hpp>
#include <redma/array/BlockVector.hpp>
#include <redma/array/BlockMatrix.hpp>

namespace RedMA {

class aPressureDrop : public aFunctionProvider {

public:
    aPressureDrop() {}

    virtual shp<aMatrix> getPressureMass(const double& time,
                                         const shp<aVector>& sol) override
        {throw new Exception("'getPressureMass' method is not implemented in class aPressureDrop");};

    virtual inline void setFlowRate(const double& Q) {M_Q = Q;}

    virtual void apply0DirichletBCs(shp<aVector> vector) const override {}

    virtual void applyDirichletBCs(const double& time, shp<aVector> vector) const override {}

    void setExtrapolatedSolution(const shp<aVector>& exSol) override
    {throw new Exception("'setExtrapolatedSolution' method must still be implemented in aPressureDrop class");}

    virtual std::map<unsigned int, std::vector<shp<aVector>>> importSolution(const std::string& filename) const override
    {throw new Exception("'importSolution' method must still be implemented in aPressureDrop class");};

protected:
    double                            M_Q;  // outflow rate

};

} // namespace RedMA


#endif //APRESSUREDROP_HPP
