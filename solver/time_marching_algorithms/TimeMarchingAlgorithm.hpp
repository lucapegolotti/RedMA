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

#include <sstream>

#include <GlobalBlockMatrix.hpp>
#include <GlobalAssembler.hpp>
#include <PrintLog.hpp>
#include <GlobalSolverOperator.hpp>
#include <GlobalSolverPreconditionerOperator.hpp>
#include <GlobalSIMPLEOperator.hpp>
#include <GlobalSIMPLEOperatorPseudoFSI.hpp>
#include <GlobalSolverOperator.hpp>
// #include <GlobalIdentityOperator.hpp>

#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/linear_algebra/InvertibleOperator.hpp>
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
    typedef Teuchos::ParameterList                          ParameterList;
    typedef std::shared_ptr<ParameterList>                  ParameterListPtr;

public:
    TimeMarchingAlgorithm(const GetPot& datafile,
                          GlobalAssemblerType* assembler,
                          commPtr_Type comm,
                          bool verbose = false);

    virtual void solveTimestep(const double &time, double &dt) = 0;

    void solveLinearSystem(GlobalBlockMatrix matrix, VectorPtr rhs,
                           VectorPtr sol);

    void solveNonLinearSystem(std::function<VectorPtr(VectorPtr)> fun,
                              std::function<GlobalBlockMatrix(VectorPtr)> jac,
                              VectorPtr sol, const double& tol,
                              const unsigned int& itMax);

    void setSolversOptions();

    void buildPreconditioner(GlobalBlockMatrix matrix);

    VectorPtr getSolution();

    unsigned int getOrder();

    void setInitialCondition(VectorPtr initalCondition);

protected:
    GetPot                                                                M_datafile;
    VectorPtr                                                             M_solution;
    GlobalAssemblerType*                                                  M_globalAssembler;
    commPtr_Type                                                          M_comm;
    bool                                                                  M_verbose;
    std::shared_ptr<LifeV::Operators::GlobalSolverOperator>               M_oper;
    std::shared_ptr<LifeV::Operators::GlobalSolverPreconditionerOperator> M_prec;
    std::shared_ptr<LifeV::Operators::InvertibleOperator>                 M_invOper;
    ParameterListPtr                                                      M_pListLinSolver;
    Teuchos::RCP<Teuchos::ParameterList>                                  M_solversOptions;
    unsigned int                                                          M_order;
};

}  // namespace RedMA

#include <TimeMarchingAlgorithm_imp.hpp>

#endif  // TIMEMARCHINGALGORITHM_HPP
