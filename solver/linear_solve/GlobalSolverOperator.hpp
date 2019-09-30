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
#define GLOBALSOLVEROPERATOR_HPP

#include <Epetra_Import.h>
#include <boost/numeric/ublas/matrix.hpp>

#include <GlobalBlockMatrix.hpp>
#include <lifev/core/linear_algebra/BlockEpetra_Map.hpp>
#include <lifev/core/linear_algebra/BlockEpetra_MultiVector.hpp>
#include <lifev/core/linear_algebra/LinearOperatorAlgebra.hpp>

namespace LifeV
{

namespace Operators
{

class GlobalSolverOperator : public LinearOperatorAlgebra
{
public:

    typedef LinearOperatorAlgebra super;
    typedef super::comm_Type comm_Type;
    typedef super::commPtr_Type commPtr_Type;
    typedef super::map_Type map_Type;
    typedef super::mapPtr_Type mapPtr_Type;
    typedef super::operator_Type operator_Type;
    typedef super::operatorPtr_Type operatorPtr_Type;
    typedef super::vector_Type vector_Type;
    typedef super::vectorPtr_Type vectorPtr_Type;

    typedef boost::numeric::ublas::matrix<operatorPtr_Type> operatorPtrContainer_Type;
    typedef std::vector<vectorPtr_Type> vectorPtrContainer_Type;
    typedef std::vector<mapPtr_Type > mapPtrContainer_Type;

    GlobalSolverOperator();

    void setUp(const std::shared_ptr<BlockEpetra_Map> & map,
               const commPtr_Type & comm);

    void setUp(const std::shared_ptr<BlockEpetra_Map> & domainMap,
               const std::shared_ptr<BlockEpetra_Map> & rangeMap,
               const commPtr_Type & comm);

    void setUp(operatorPtrContainer_Type blockOper,
               const commPtr_Type & comm);

    void setBlock(UInt iblock, UInt jblock, const operatorPtr_Type & operBlock);

    void fillComplete();

    int SetUseTranspose(bool useTranspose);

    virtual int Apply(const vector_Type & X, vector_Type & Y) const;

    virtual int ApplyInverse(const vector_Type & X, vector_Type & Y) const;

    double NormInf() const {return -1;}

    virtual const char * Label() const {return M_name.c_str();}

    bool UseTranspose() const {return M_useTranspose;}

    bool HasNormInf() const {return false;}

    const comm_Type & Comm() const {return *M_comm;}

    const operatorPtr_Type& block (UInt iblock, UInt jblock) const;

    const map_Type & OperatorDomainMap() const {return *(M_domainMap->monolithicMap());}

    const mapPtr_Type & OperatorDomainMap_ptr() const {return M_domainMap->monolithicMap();}

    const std::shared_ptr<BlockEpetra_Map> & OperatorDomainBlockMapPtr() const {return M_domainMap;}

    const map_Type & OperatorRangeMap() const {return *(M_rangeMap->monolithicMap());}

    const mapPtr_Type & OperatorRangeMap_ptr() const {return M_rangeMap->monolithicMap();}

    const std::shared_ptr<BlockEpetra_Map> & OperatorRangeBlockMapPtr() const {return M_rangeMap;}

protected:
    int applyNoTranspose(const vector_Type & X, vector_Type & Y) const;

    int applyTranspose(const vector_Type & X, vector_Type & Y) const;

    int blockJacobi(const vector_Type & X, vector_Type & Y) const;

    int blockUpperTriangularSolve(const vector_Type & X, vector_Type & Y) const;

    int blockLowerTriangularSolve(const vector_Type & X, vector_Type & Y) const;

    void setName(const std::string & name){M_name =name;}

private:

    UInt M_nBlockRows;

    UInt M_nBlockCols;

    std::string M_name;

    commPtr_Type M_comm;

    std::shared_ptr<BlockEpetra_Map> M_domainMap;

    std::shared_ptr<BlockEpetra_Map> M_rangeMap;

    operatorPtrContainer_Type M_oper;

    bool M_useTranspose;
};

}

}
#endif // GLOBALSOLVEROPERATOR_HPP
