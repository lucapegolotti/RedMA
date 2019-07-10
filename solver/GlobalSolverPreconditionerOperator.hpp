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

#ifndef GLOBALSOLVERPRECONDITIONEROPERATOR_HPP
#define GLOBALSOLVERPRECONDITIONEROPERATOR_HPP 1

#include <boost/numeric/ublas/matrix.hpp>

#include <lifev/core/linear_algebra/LinearOperatorAlgebra.hpp>
#include <lifev/core/util/Factory.hpp>
#include <lifev/core/util/FactorySingleton.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/linear_algebra/BlockEpetra_Map.hpp>

#include <GlobalBlockMatrix.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV
{
namespace Operators
{

class GlobalSolverPreconditionerOperator : public LinearOperatorAlgebra
{
public:

    typedef MatrixEpetra<Real> matrixEpetra_Type;
    typedef std::shared_ptr<matrixEpetra_Type> matrixEpetraPtr_Type;
    typedef Epetra_Comm comm_Type;
    typedef std::shared_ptr<comm_Type> commPtr_Type;
    typedef Epetra_Map map_Type;
    typedef std::shared_ptr<map_Type> mapPtr_Type;
    typedef std::shared_ptr<const map_Type> constMapPtr_Type;
    typedef Epetra_Operator operator_Type;
    typedef std::shared_ptr<operator_Type> operatorPtr_Type;
    typedef  std::shared_ptr<Teuchos::ParameterList> parameterListPtr_Type;

    typedef boost::numeric::ublas::matrix<operatorPtr_Type> operatorPtrContainer_Type;
    typedef std::vector<vectorPtr_Type>                     vectorPtrContainer_Type;
    typedef std::vector<mapPtr_Type >                       mapPtrContainer_Type;

    GlobalSolverPreconditionerOperator();

    ~GlobalSolverPreconditionerOperator();

    virtual int SetUseTranspose(bool UseTranspose){};

    virtual int Apply(const vector_Type& X, vector_Type& Y) const {};

    virtual int ApplyInverse(const vector_Type& X, vector_Type& Y) const {};

    virtual void setUp(RedMA::GlobalBlockMatrix matrix, const commPtr_Type& comm) = 0;

    void setSolversOptions(const Teuchos::ParameterList& solversOptions)
    {
        M_solversOptions = solversOptions;
    }

    virtual double NormInf() const {};

    virtual const char * Label() const {};

    virtual bool UseTranspose() const {};

    virtual bool HasNormInf() const {};

    virtual const comm_Type & Comm() const {};

    virtual const map_Type & OperatorDomainMap() const {};

    virtual const map_Type & OperatorRangeMap() const {};

protected:
    Teuchos::ParameterList M_solversOptions;
};
inline GlobalSolverPreconditionerOperator::GlobalSolverPreconditionerOperator()
{
}
inline GlobalSolverPreconditionerOperator::~GlobalSolverPreconditionerOperator()
{
}
typedef FactorySingleton<Factory<GlobalSolverPreconditionerOperator, std::string> >
                                                GlobalPreconditionerFactory;
}
}
#endif // GLOBALSOLVERPRECONDITIONEROPERATOR_HPP
