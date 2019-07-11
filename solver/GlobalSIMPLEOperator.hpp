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

#ifndef GLOBALSIMPLEOPERATOR_H
#define GLOBALSIMPLEOPERATOR_H 1

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>

#include <lifev/core/linear_algebra/BlockEpetra_Map.hpp>

#include <boost/numeric/ublas/matrix.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <lifev/core/linear_algebra/LinearOperatorAlgebra.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>

#include <lifev/navier_stokes_blocks/solver/aSIMPLEOperator.hpp>

#include <GlobalSolverPreconditionerOperator.hpp>

namespace LifeV
{
namespace Operators
{

class GlobalSIMPLEOperator : public GlobalSolverPreconditionerOperator
{
public:

    typedef LinearOperatorAlgebra                           super;
    typedef super::comm_Type                                comm_Type;
    typedef super::commPtr_Type                             commPtr_Type;
    typedef super::map_Type                                 map_Type;
    typedef super::mapPtr_Type                              mapPtr_Type;
    typedef super::operator_Type                            operator_Type;
    typedef super::operatorPtr_Type                         operatorPtr_Type;
    typedef super::vector_Type                              vector_Type;
    typedef super::vectorPtr_Type                           vectorPtr_Type;
    typedef MatrixEpetra<Real>                              matrixEpetra_Type;
    typedef std::shared_ptr<matrixEpetra_Type>              matrixEpetraPtr_Type;
    typedef VectorEpetra                                    vectorEpetra_Type;
    typedef std::shared_ptr<vectorEpetra_Type>              vectorEpetraPtr_Type;
    typedef LifeV::MapEpetra                                mapEpetra_Type;
    typedef std::shared_ptr<mapEpetra_Type>                 mapEpetraPtr_Type;

    typedef boost::numeric::ublas::matrix<operatorPtr_Type> operatorPtrContainer_Type;
    typedef std::shared_ptr<Operators::ApproximatedInvertibleRowMatrix>
                                                  ApproximatedInvertedMatrixPtr;
    typedef boost::numeric::ublas::matrix<ApproximatedInvertedMatrixPtr>
                                           GridApproximatedInvertedMatricesPtrs;
    typedef std::vector<vectorPtr_Type>                     vectorPtrContainer_Type;
    typedef std::vector<mapPtr_Type >                       mapPtrContainer_Type;

    typedef std::shared_ptr<LifeV::Operators::NavierStokesPreconditionerOperator>
                                                            PreconditionerPtr;

    GlobalSIMPLEOperator();

    virtual ~GlobalSIMPLEOperator();

    void setUp(RedMA::GlobalBlockMatrix matrix, const commPtr_Type & comm);

    int SetUseTranspose(bool UseTranspose){M_useTranspose = UseTranspose; return 0;}

    void computeBAm1BT_inverse(unsigned int rowIndex, unsigned int colIndex);

    int Apply(const vector_Type &/*X*/, vector_Type &/*Y*/) const {return -1;};

    int ApplyInverse(const vector_Type &X, vector_Type &Y) const;

    double NormInf() const {return -1.0;}

    void fillComplete();

    const char * Label() const {return M_label.c_str();}

    bool UseTranspose() const {return M_useTranspose;}

    bool HasNormInf() const {return false;}

    const comm_Type & Comm() const {return *M_comm;}

    void computeAllBAm1Binverses();

private:

    std::vector<PreconditionerPtr>       M_SIMPLEOperators;
    RedMA::GlobalBlockMatrix             M_matrix;
    UInt                                 M_nBlockRows;
    UInt                                 M_nBlockCols;
    std::string                          M_name;
    commPtr_Type                         M_comm;
    std::shared_ptr<BlockEpetra_Map>     M_domainMap;
    std::shared_ptr<BlockEpetra_Map>     M_rangeMap;
    operatorPtrContainer_Type            M_oper;
    std::string                          M_label;
    bool                                 M_useTranspose;
    unsigned int                         M_nPrimalBlocks;
    GridApproximatedInvertedMatricesPtrs M_approximatedBAm1Binverses;
};

inline GlobalSolverPreconditionerOperator * create_GlobalSIMPLE()
{
    return new GlobalSIMPLEOperator ();
}
namespace
{
static bool S_register_aGlobalSimple =
GlobalPreconditionerFactory::instance().registerProduct("GlobalSIMPLE", &create_GlobalSIMPLE);
}
}
}
#endif // GLOBALSIMPLEOPERATOR_H
