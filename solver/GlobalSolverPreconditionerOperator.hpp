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

#ifndef GLOBALSOLVEROPERATOR_HPP
#define GLOBALSOLVEROPERATOR_HPP 1

#include <lifev/core/linear_algebra/LinearOperatorAlgebra.hpp>
#include <lifev/core/util/Factory.hpp>
#include <lifev/core/util/FactorySingleton.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/linear_algebra/BlockEpetra_Map.hpp>

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

    GlobalSolverPreconditionerOperator();

    ~GlobalSolverPreconditionerOperator();

    virtual int SetUseTranspose(bool UseTranspose){};

    virtual void SetSolveMomentumBelos(){};

    virtual int Apply(const vector_Type& X, vector_Type& Y) const {};


    virtual int ApplyInverse(const vector_Type& X, vector_Type& Y) const {};

    virtual double NormInf() const {};

    virtual const char * Label() const {};

    virtual bool UseTranspose() const {};

    virtual bool HasNormInf() const {};

    virtual const comm_Type & Comm() const {};

    virtual const map_Type & OperatorDomainMap() const {};

    virtual const map_Type & OperatorRangeMap() const {};

    virtual void setUp (const matrixEpetraPtr_Type & F,
       	   	            const matrixEpetraPtr_Type & B,
       	   	            const matrixEpetraPtr_Type & Btranspose){};

    virtual void setUp (const matrixEpetraPtr_Type & F,
                        const matrixEpetraPtr_Type & B,
                        const matrixEpetraPtr_Type & Btranspose,
                        const matrixEpetraPtr_Type & D){};

    virtual void setUp(const matrixEpetraPtr_Type & F,
                       const matrixEpetraPtr_Type & B,
                       const matrixEpetraPtr_Type & Btranspose,
                       const matrixEpetraPtr_Type & Fp,
                       const matrixEpetraPtr_Type & Mp,
                       const matrixEpetraPtr_Type & Mu){};

    virtual void setVelocityMass(const matrixEpetraPtr_Type & Mu){ };

    virtual void setPressureMass(const matrixEpetraPtr_Type & Mp) { };

    virtual void setOptions(const Teuchos::ParameterList& solversOptions){};

    virtual void setDomainMap(const std::shared_ptr<BlockEpetra_Map> & domainMap){};

    virtual void setRangeMap(const std::shared_ptr<BlockEpetra_Map> & rangeMap){};

    virtual void updateApproximatedMomentumOperator(){};

    virtual void updateApproximatedSchurComplementOperator(){};

    virtual void updateApproximatedPressureMassOperator(){};

    virtual void setMomentumOptions(const parameterListPtr_Type & _oList){};

    virtual void setSchurOptions(const parameterListPtr_Type & _oList){};

    virtual void setPressureMassOptions(const parameterListPtr_Type & _oList){};

    virtual matrixEpetraPtr_Type const& F() const {};

    virtual matrixEpetraPtr_Type const& B() const {};

    virtual matrixEpetraPtr_Type const& Btranspose() const {};

    virtual int ApplyInverse( VectorEpetra const& X_velocity,
    						  VectorEpetra const& X_pressure,
    						  VectorEpetra & Y_velocity,
    						  VectorEpetra & Y_pressure) const {};

    virtual int setCollectSnapshot( bool _collectRbSnapshot ) { };

    virtual int setUseRbApproximation( bool _useRbApproximation ) { };

    virtual int setParam( const Epetra_SerialDenseVector& _mu ){};

    virtual double getAvgVelocityIteration( ) { };

    virtual int resetIterations( ) { };

    virtual int setViscosity( Real viscosity ) { };

private:

};
inline GlobalSolverPreconditionerOperator::GlobalSolverPreconditionerOperator()
{
}
inline GlobalSolverPreconditionerOperator::~GlobalSolverPreconditionerOperator()
{
}
typedef FactorySingleton<Factory<GlobalSolverPreconditionerOperator, std::string> >
                                                        NSPreconditionerFactory;
}
}
#endif // GLOBALSOLVEROPERATOR_HPP
