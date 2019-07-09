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

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>

#include <GlobalSolverPreconditionerOperator.hpp>

#ifndef GLOBALIDENTITYOPERATOR_H
#define GLOBALIDENTITYOPERATOR_H 1

namespace LifeV
{
namespace Operators
{

class GlobalIdentityOperator: public GlobalSolverPreconditionerOperator
{
public:

    typedef  Epetra_MultiVector                        vector_Type;
    typedef  std::shared_ptr<vector_Type>              vectorPtr_Type;
    typedef  Epetra_Map                                map_Type;
    typedef  std::shared_ptr<map_Type>                 mapPtr_Type;
    typedef  LinearOperatorAlgebra                     super;
    typedef  Epetra_CrsMatrix                          matrix_Type;
    typedef  std::shared_ptr<matrix_Type>              matrixPtr_Type;
    typedef  MatrixEpetra<Real>                        matrixEpetra_Type;
    typedef  std::shared_ptr<matrixEpetra_Type>        matrixEpetraPtr_Type;
    typedef  Epetra_Vector                             lumpedMatrix_Type;
    typedef  std::shared_ptr<lumpedMatrix_Type>        lumpedMatrixPtr_Type;
    typedef  super::comm_Type                          comm_Type;
    typedef  super::commPtr_Type                       commPtr_Type;
    typedef  std::shared_ptr<Teuchos::ParameterList>   parameterListPtr_Type;
    typedef  MapEpetra                                 mapEpetra_Type;
    typedef  std::shared_ptr<mapEpetra_Type>           mapEpetraPtr_Type;
    typedef  VectorEpetra                              VectorEpetra_Type;
    typedef  std::shared_ptr<VectorEpetra_Type>        VectorEpetraPtr_Type;
    typedef  LifeV::Preconditioner                     prec_Type;
    typedef  std::shared_ptr<prec_Type>                precPtr_Type;

    GlobalIdentityOperator();

    virtual ~GlobalIdentityOperator();

    virtual void setUp(operatorPtrContainer_Type oper, const commPtr_Type& comm);

    int SetUseTranspose(bool UseTranspose){M_useTranspose = UseTranspose; return 0;}

    void setDomainMap(const std::shared_ptr<BlockEpetra_Map> & domainMap)
    {
        M_operatorDomainMap = domainMap;
    }

    void setRangeMap(const std::shared_ptr<BlockEpetra_Map> & rangeMap)
    {
        M_operatorRangeMap = rangeMap;
    }

    int Apply(const vector_Type &/*X*/, vector_Type &/*Y*/) const {return -1;};

    int ApplyInverse(const vector_Type &X, vector_Type &Y) const;

    double NormInf() const {return -1.0;}

    const char * Label() const {return M_label.c_str();}

    bool UseTranspose() const {return M_useTranspose;}

    bool HasNormInf() const {return false;}

    const comm_Type & Comm() const {return *M_comm;}

    const map_Type & OperatorDomainMap() const
    {
        return *(M_operatorDomainMap->monolithicMap());
    }

    const map_Type & OperatorRangeMap() const
    {
        return *(M_operatorRangeMap->monolithicMap());
    }

    void showMe();

private:

    std::shared_ptr<BlockEpetra_Map> M_operatorDomainMap;

    std::shared_ptr<BlockEpetra_Map> M_operatorRangeMap;\

    commPtr_Type M_comm;

    bool M_useTranspose;

    const std::string M_label;
};

inline GlobalSolverPreconditionerOperator * create_Identity()
{
    return new GlobalIdentityOperator ();
}
namespace
{
static bool S_register_aSimple =
GlobalPreconditionerFactory::instance().registerProduct("SIMPLE", &create_Identity);
}
}
}
#endif // GLOBALIDENTITYOPERATOR_H
