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
    using TimeMarchingAlgorithm<AssemblerType>::M_datafile;
    using TimeMarchingAlgorithm<AssemblerType>::M_solution;
    using TimeMarchingAlgorithm<AssemblerType>::M_globalAssembler;
    using TimeMarchingAlgorithm<AssemblerType>::solveLinearSystem;
public:
    typedef typename TimeMarchingAlgorithm<AssemblerType>::GlobalAssemblerType
                     GlobalAssemblerType;
    typedef typename TimeMarchingAlgorithm<AssemblerType>::Vector
                     Vector;
    typedef typename TimeMarchingAlgorithm<AssemblerType>::VectorPtr
                     VectorPtr;
    typedef typename TimeMarchingAlgorithm<AssemblerType>::Matrix
                     Matrix;
    typedef typename TimeMarchingAlgorithm<AssemblerType>::MatrixPtr
                     MatrixPtr;
    typedef typename TimeMarchingAlgorithm<AssemblerType>::MapEpetra
                     MapEpetra;
    typedef typename TimeMarchingAlgorithm<AssemblerType>::MapEpetraPtr
                     MapEpetraPtr;
    typedef typename TimeMarchingAlgorithm<AssemblerType>::commPtr_Type
                     commPtr_Type;

    RosenbrockAlgorithm(const GetPot& datafile,
                        GlobalAssemblerType* assembler,
                        commPtr_Type comm);

    virtual void solveTimestep(const double &time, double &dt);

private:
    LifeV::RosenbrockCoeff  M_coefficients;
};

}  // namespace RedMA

#include <RosenbrockAlgorithm_imp.hpp>

#endif  // ROSENBROCKALGORITHM_HPP
