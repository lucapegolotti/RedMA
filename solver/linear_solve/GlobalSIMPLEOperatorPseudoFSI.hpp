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

#ifndef GLOBALSIMPLEOPERATORPSEUDOFSI_H
#define GLOBALSIMPLEOPERATORPSEUDOFSI_H 1

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>

#include <lifev/core/linear_algebra/BlockEpetra_Map.hpp>

#include <boost/numeric/ublas/matrix.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <lifev/core/linear_algebra/LinearOperatorAlgebra.hpp>

#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>
#include <lifev/core/array/MatrixEpetraStructuredUtility.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>

#include <lifev/core/util/LifeChrono.hpp>

#include <PrintLog.hpp>

#include <lifev/navier_stokes_blocks/solver/aSIMPLEOperator.hpp>

#include <GlobalSolverPreconditionerOperator.hpp>

namespace LifeV
{
namespace Operators
{

class GlobalSIMPLEOperatorPseudoFSI : public GlobalSolverPreconditionerOperator
{
public:

    typedef LinearOperatorAlgebra                       super;
    typedef super::comm_Type                            comm_Type;
    typedef super::commPtr_Type                         commPtr_Type;
    typedef super::map_Type                             map_Type;
    typedef super::mapPtr_Type                          mapPtr_Type;
    typedef super::operator_Type                        operator_Type;
    typedef super::operatorPtr_Type                     operatorPtr_Type;
    typedef super::vector_Type                          vector_Type;
    typedef super::vectorPtr_Type                       vectorPtr_Type;
    typedef MatrixEpetra<Real>                          matrixEpetra_Type;
    typedef std::shared_ptr<matrixEpetra_Type>          matrixEpetraPtr_Type;
    typedef VectorEpetra                                vectorEpetra_Type;
    typedef std::shared_ptr<vectorEpetra_Type>          vectorEpetraPtr_Type;
    typedef LifeV::MapEpetra                            mapEpetra_Type;
    typedef std::shared_ptr<mapEpetra_Type>             mapEpetraPtr_Type;

    typedef boost::numeric::ublas::matrix<operatorPtr_Type>
                                                  operatorPtrContainer_Type;
    typedef Operators::ApproximatedInvertibleRowMatrix
                                                  ApproximatedInvertibleMatrix;
    typedef std::shared_ptr<ApproximatedInvertibleMatrix>
                                                  ApproximatedInvertibleMatrixPtr;
    typedef std::vector<vectorPtr_Type>                 vectorPtrContainer_Type;
    typedef std::vector<mapPtr_Type >                   mapPtrContainer_Type;

    typedef std::shared_ptr<LifeV::Operators::NavierStokesPreconditionerOperator>
                                                            PreconditionerPtr;

    GlobalSIMPLEOperatorPseudoFSI();

    virtual ~GlobalSIMPLEOperatorPseudoFSI();

    void setUp(RedMA::BlockMatrix matrix, const commPtr_Type & comm);

    int SetUseTranspose(bool UseTranspose){M_useTranspose = UseTranspose; return 0;}

    void computeAm1BT(unsigned int rowIndex, unsigned int colIndex);

    int Apply(const vector_Type &/*X*/, vector_Type &/*Y*/) const {return -1;};

    int ApplyInverse(const vector_Type &X, vector_Type &Y) const;

    double NormInf() const {return -1.0;}

    void fillComplete();

    const char * Label() const {return M_label.c_str();}

    bool UseTranspose() const {return M_useTranspose;}

    bool HasNormInf() const {return false;}

    const comm_Type & Comm() const {return *M_comm;}

    void computeGlobalSchurComplement();

    void applyEverySIMPLEOperator(const vectorEpetra_Type& X,
                                  vectorEpetra_Type &Y) const;

    void applyEveryB(const vectorEpetra_Type& X, vectorEpetra_Type &Y) const;

    void applyEveryBT(const vectorEpetra_Type& X, vectorEpetra_Type &Y) const;

private:

    std::vector<PreconditionerPtr>       M_SIMPLEOperators;
    RedMA::BlockMatrix             M_matrix;
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
    RedMA::BlockMatrix             M_Am1BT;
    mapEpetraPtr_Type                    M_monolithicMap;
    std::vector<mapEpetraPtr_Type>       M_allMaps;
    mapEpetraPtr_Type                    M_primalMap;
    std::vector<mapEpetraPtr_Type>       M_primalMaps;
    mapEpetraPtr_Type                    M_dualMap;
    std::vector<mapEpetraPtr_Type>       M_dualMaps;
    matrixEpetraPtr_Type                 M_globalSchurComplement;
    std::vector<unsigned int>            M_dimensionsInterfaces;
    ApproximatedInvertibleMatrixPtr      M_approximatedGlobalSchurInverse;
    std::vector<double>                  M_blockNorm22;
    std::vector<double>                  M_blockNorm20;
};

inline GlobalSolverPreconditionerOperator * create_GlobalSIMPLEPseudoFSI()
{
    return new GlobalSIMPLEOperatorPseudoFSI ();
}
namespace
{
static bool S_register_aGlobalSimplePseudoFSI =
GlobalPreconditionerFactory::instance().registerProduct("GlobalSIMPLEPseudoFSI", &create_GlobalSIMPLEPseudoFSI);
}
}
}
#endif // GLOBALSIMPLEOPERATORPSEUDOFSI_H
