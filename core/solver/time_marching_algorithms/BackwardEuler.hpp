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

#ifndef BACKWARDEULER_HPP
#define BACKWARDEULER_HPP

#include <TimeMarchingAlgorithm.hpp>
#include <PrintLog.hpp>
#include <lifev/core/array/VectorBlockStructure.hpp>
#include <lifev/core/array/MatrixBlockStructure.hpp>
#include <BlockMatrix.hpp>

#include <lifev/navier_stokes_blocks/solver/RosenbrockCoeff.hpp>

namespace RedMA
{

class BackwardEuler : public TimeMarchingAlgorithm
{
    using TimeMarchingAlgorithm::M_datafile;
    using TimeMarchingAlgorithm::M_solution;
    using TimeMarchingAlgorithm::M_globalAssembler;
    using TimeMarchingAlgorithm::M_verbose;
    using TimeMarchingAlgorithm::solveLinearSystem;
    using TimeMarchingAlgorithm::getSolution;
public:
    typedef typename TimeMarchingAlgorithm::GlobalAssemblerType
                     GlobalAssemblerType;
    typedef typename TimeMarchingAlgorithm::Vector
                     Vector;
    typedef typename TimeMarchingAlgorithm::VectorPtr
                     VectorPtr;
    typedef typename TimeMarchingAlgorithm::Matrix
                     Matrix;
    typedef typename TimeMarchingAlgorithm::MatrixPtr
                     MatrixPtr;
    typedef typename TimeMarchingAlgorithm::MapEpetra
                     MapEpetra;
    typedef typename TimeMarchingAlgorithm::MapEpetraPtr
                     MapEpetraPtr;
    typedef typename TimeMarchingAlgorithm::commPtr_Type
                     commPtr_Type;

    BackwardEuler(const GetPot& datafile,
                        GlobalAssemblerType* assembler,
                        commPtr_Type comm,
                        bool verbose = false);

    virtual void solveTimestep(const double &time, double &dt);

    VectorPtr assembleF(const double& time, VectorPtr tentativeSol,
                        const double& dt);

    BlockMatrix assembleJac(const double& time, VectorPtr tentativeSol,
                                  const double& dt);

private:
    VectorPtr               M_prevSolution;
};

}  // namespace RedMA

#endif  // BACKWARDEULER_HPP
