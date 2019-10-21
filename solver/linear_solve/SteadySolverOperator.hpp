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

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>

#include <lifev/core/linear_algebra/BlockEpetra_Map.hpp>
#include <lifev/navier_stokes_blocks/solver/NavierStokesPreconditionerOperator.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>

#ifndef STEADYSOLVEROPERATOR_HPP
#define STEADYSOLVEROPERATOR_HPP 1

namespace LifeV
{
namespace Operators
{

class SteadySolverOperator : public NavierStokesPreconditionerOperator
{
public:
    //! @name Public Types
    //@{

    typedef  Epetra_MultiVector                        vector_Type;
    typedef  std::shared_ptr<vector_Type>            vectorPtr_Type;
    typedef  Epetra_Map                                map_Type;
    typedef  std::shared_ptr<map_Type> 			   mapPtr_Type;
    typedef  LinearOperatorAlgebra                     super;
    typedef  Epetra_CrsMatrix                          matrix_Type;
    typedef  std::shared_ptr<matrix_Type>            matrixPtr_Type;
    typedef  MatrixEpetra<Real>                        matrixEpetra_Type;
    typedef  std::shared_ptr<matrixEpetra_Type>      matrixEpetraPtr_Type;
    typedef  Epetra_Vector                             lumpedMatrix_Type;
    typedef  std::shared_ptr<lumpedMatrix_Type>      lumpedMatrixPtr_Type;
    typedef  super::comm_Type                          comm_Type;
    typedef  super::commPtr_Type                       commPtr_Type;
    typedef  std::shared_ptr<Teuchos::ParameterList> parameterListPtr_Type;
    typedef  MapEpetra                                 mapEpetra_Type;
    typedef  std::shared_ptr<mapEpetra_Type>         mapEpetraPtr_Type;
    typedef  VectorEpetra                              VectorEpetra_Type;
    typedef  std::shared_ptr<VectorEpetra_Type>      VectorEpetraPtr_Type;
    typedef  LifeV::Preconditioner                     prec_Type;
    typedef  std::shared_ptr<prec_Type>              precPtr_Type;

    SteadySolverOperator();

    virtual ~SteadySolverOperator();

    void setUp(const matrixEpetraPtr_Type & F,
               const matrixEpetraPtr_Type & B,
               const matrixEpetraPtr_Type & Btranspose);

    void setUp (const matrixEpetraPtr_Type & F,
                const matrixEpetraPtr_Type & B,
                const matrixEpetraPtr_Type & M_p,
                const matrixEpetraPtr_Type & D);

    int SetUseTranspose(bool UseTranspose){M_useTranspose = UseTranspose; return 0;}
    //! set the domain map
    void setDomainMap(const std::shared_ptr<BlockEpetra_Map> & domainMap){M_operatorDomainMap = domainMap;}
    //! set the range map
    void setRangeMap(const std::shared_ptr<BlockEpetra_Map> & rangeMap){M_operatorRangeMap = rangeMap;}

    int Apply(const vector_Type &/*X*/, vector_Type &/*Y*/) const {return -1;};

    int ApplyInverse( VectorEpetra_Type const& X_velocity,
                      VectorEpetra_Type const& X_pressure,
                      VectorEpetra_Type & Y_velocity,
                      VectorEpetra_Type & Y_pressure) const;

    int ApplyInverse(const vector_Type &X, vector_Type &Y) const;

    double NormInf() const {return -1.0;}

    const char * Label() const {return M_label.c_str();}

    bool UseTranspose() const {return M_useTranspose;}

    bool HasNormInf() const {return false;}

    const comm_Type & Comm() const {return *M_comm;}

    const map_Type & OperatorDomainMap() const {return *(M_operatorDomainMap->monolithicMap());}

    const map_Type & OperatorRangeMap() const {return *(M_operatorRangeMap->monolithicMap());}

    void updateApproximatedMomentumOperator();

    void updateApproximatedPressureMassOperator();

    void setMomentumOptions(const parameterListPtr_Type & _oList);

    void setPressureMassOptions(const parameterListPtr_Type & _oList);

    void SetSolveMomentumBelos(){M_solveMomentumBelos = true;};

    void showMe();

    void setOptions(const Teuchos::ParameterList& solversOptions);

    matrixEpetraPtr_Type const& F() const { return M_F; }

    matrixEpetraPtr_Type const& Btranspose() const { return M_Btranspose; }

private:
    void setMaps();

    std::shared_ptr<BlockEpetra_Map> M_operatorDomainMap;

    std::shared_ptr<BlockEpetra_Map> M_operatorRangeMap;

    matrixEpetraPtr_Type M_F;

    matrixEpetraPtr_Type M_Btranspose;

    matrixEpetraPtr_Type M_Mp;

    commPtr_Type M_comm;

    bool M_useTranspose;

    std::shared_ptr<Operators::ApproximatedInvertibleRowMatrix> M_approximatedMomentumOperator;

    std::shared_ptr<Operators::ApproximatedInvertibleRowMatrix> M_approximatedPressureOperator;

    parameterListPtr_Type M_momentumOptions;

    parameterListPtr_Type M_pressureOptions;

    mapEpetraPtr_Type M_monolithicMap;

    bool M_solveMomentumBelos;

    //! Label
    const std::string M_label;

    //! Vectors needed for the apply inverse
    std::shared_ptr<VectorEpetra_Type> M_Z;

    std::shared_ptr<Epetra_Vector> M_invD;

    std::shared_ptr<VectorEpetra_Type> M_X_velocity;
    std::shared_ptr<VectorEpetra_Type> M_X_pressure;
    std::shared_ptr<VectorEpetra_Type> M_Y_velocity;
    std::shared_ptr<VectorEpetra_Type> M_Y_pressure;

    bool M_useStabilization;

    std::shared_ptr<LinearSolver>               M_linearSolver;
    precPtr_Type                                M_vStiffnessPreconditioner;
    Teuchos::RCP< Teuchos::ParameterList >      M_belosList;
};

//! Factory create function
inline NavierStokesPreconditionerOperator * create_SteadySolver()
{
    return new SteadySolverOperator();
}
namespace
{
static bool S_register_SteadySolver = NSPreconditionerFactory::instance().registerProduct("STEADY", &create_SteadySolver);
}

} /* end namespace Operators */
} //end namespace
#endif // STEADYSOLVEROPERATOR_HPP
