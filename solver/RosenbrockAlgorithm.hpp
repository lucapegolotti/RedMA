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

#ifndef ROSENBROCKALGORITHM_HPP
#define ROSENBROCKALGORITHM_HPP

#include <TimeMarchingAlgorithm.hpp>
#include <LinearSolver.hpp>

#include <lifev/core/array/VectorBlockStructure.hpp>
#include <lifev/core/array/MatrixBlockStructure.hpp>

#include <lifev/navier_stokes_blocks/solver/RosenbrockCoeff.hpp>

namespace RedMA
{

template <class AssemblerType>
class RosenbrockAlgorithm : public TimeMarchingAlgorithm<AssemblerType>
{
protected:
    typedef GlobalAssembler<AssemblerType>              GlobalAssemblerType;
    typedef LifeV::VectorEpetra                         Vector;
    typedef std::shared_ptr<Vector>                     VectorPtr;
    typedef LifeV::MatrixEpetra<double>                 Matrix;
    typedef std::shared_ptr<Matrix>                     MatrixPtr;
    typedef LifeV::MapEpetra                            MapEpetra;
    typedef std::shared_ptr<MapEpetra>                  MapEpetraPtr;

public:
    RosenbrockAlgorithm(const GetPot& datafile);

    virtual void solveTimestep(const double &time, double &dt,
                               GlobalAssemblerType& assembler,
                               const LinearSolver& linearSolver);

private:
    LifeV::RosenbrockCoeff  M_coefficients;
    VectorPtr               M_solution;
};

}  // namespace RedMA

#include <RosenbrockAlgorithm_imp.hpp>

#endif  // ROSENBROCKALGORITHM_HPP
