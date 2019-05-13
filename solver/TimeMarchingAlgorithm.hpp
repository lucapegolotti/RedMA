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

#ifndef TIMEMARCHINGALGORITHM_HPP
#define TIMEMARCHINGALGORITHM_HPP

#include <GlobalAssembler.hpp>
#include <PrintLog.hpp>

#include <lifev/core/algorithm/SolverAztecOO.hpp>

namespace RedMA
{

template <class AssemblerType>
class TimeMarchingAlgorithm
{
protected:
    typedef GlobalAssembler<AssemblerType>                  GlobalAssemblerType;
    typedef LifeV::VectorEpetra                             Vector;
    typedef std::shared_ptr<Vector>                         VectorPtr;
    typedef std::shared_ptr<Epetra_Comm>                    commPtr_Type;
    typedef LifeV::MatrixEpetra<double>                     Matrix;
    typedef std::shared_ptr<Matrix>                         MatrixPtr;
    typedef LifeV::MapEpetra                                MapEpetra;
    typedef std::shared_ptr<MapEpetra>                      MapEpetraPtr;

public:
    TimeMarchingAlgorithm(const GetPot& datafile,
                          GlobalAssemblerType* assembler,
                          commPtr_Type comm,
                          bool verbose = false);

    virtual void solveTimestep(const double &time, double &dt) = 0;

    void solveLinearSystem(MatrixPtr matrix, VectorPtr rhs, VectorPtr sol);

    void solveNonLinearSystem(std::function<VectorPtr(VectorPtr)> fun,
                              std::function<MatrixPtr(VectorPtr)> jac,
                              VectorPtr sol, const double& tol,
                              const unsigned int& itMax);

    VectorPtr getSolution();

protected:
    GetPot                  M_datafile;
    VectorPtr               M_solution;
    GlobalAssemblerType*    M_globalAssembler;
    commPtr_Type            M_comm;
    bool                    M_verbose;
};

}  // namespace RedMA

#include <TimeMarchingAlgorithm_imp.hpp>

#endif  // TIMEMARCHINGALGORITHM_HPP
