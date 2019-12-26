#include <SteadySolverOperator.hpp>

#include <lifev/core/linear_algebra/IfpackPreconditioner.hpp>
#include <lifev/core/linear_algebra/MLPreconditioner.hpp>
#include <lifev/core/linear_algebra/TwoLevelPreconditioner.hpp>
#include <lifev/core/linear_algebra/AztecooOperatorAlgebra.hpp>
#include <lifev/core/linear_algebra/BelosOperatorAlgebra.hpp>
#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>

#include <lifev/core/algorithm/PreconditionerML.hpp>

namespace LifeV
{
namespace Operators
{

SteadySolverOperator::
SteadySolverOperator():
  M_label("SteadySolverOperator"),
  M_useTranspose(false),
  M_approximatedMomentumOperator(new Operators::ApproximatedInvertibleRowMatrix),
  M_approximatedPressureOperator(new Operators::ApproximatedInvertibleRowMatrix),
  M_useStabilization(false),
  M_solveMomentumBelos(false)
{
}

SteadySolverOperator::
~SteadySolverOperator()
{
}

void
SteadySolverOperator::
showMe(){
}

void
SteadySolverOperator::setUp(const matrixEpetraPtr_Type & F,
                            const matrixEpetraPtr_Type & Btranspose,
                            const matrixEpetraPtr_Type & Mp)
{
    M_F = F;
    M_Btranspose = Btranspose;
    M_Mp = Mp;
    M_comm = F->map().commPtr();
    M_monolithicMap.reset(new mapEpetra_Type(M_F->map()));
    *M_monolithicMap += Mp->map();
    M_Z.reset(new VectorEpetra_Type(F->map(), Unique));
    M_X_velocity.reset(new VectorEpetra_Type(M_F->map(), Unique));
    M_X_pressure.reset(new VectorEpetra_Type(Mp->map(), Unique));
    M_Y_velocity.reset(new VectorEpetra_Type(M_F->map(), Unique));
    M_Y_pressure.reset(new VectorEpetra_Type(Mp->map(), Unique));
    M_useStabilization = false;

    if (M_solveMomentumBelos)
    {

        typedef LifeV::PreconditionerML                                         precML_type;
        typedef std::shared_ptr<PreconditionerML>                             precMLPtr_type;

        typedef LifeV::PreconditionerIfpack                                     precIf_type;
        typedef std::shared_ptr<precIf_type>                                  precIfPtr_type;

        precML_type * precRawPtr;
        precRawPtr = new precML_type;

        // precIf_type * precRawPtr;
        // precRawPtr = new precIf_type;

        GetPot dataFile; // fake datafile
        precRawPtr->setDataFromGetPot(dataFile, "precMLL");
        M_vStiffnessPreconditioner.reset(precRawPtr);
        M_vStiffnessPreconditioner->buildPreconditioner(M_F);

        M_belosList = Teuchos::rcp (new Teuchos::ParameterList);
        M_belosList = Teuchos::getParametersFromXmlFile("SolverParamList2.xml");
        M_linearSolver.reset(new LinearSolver(M_comm));
        M_linearSolver->setOperator(M_F);
        M_linearSolver->setParameters(*M_belosList );
        M_linearSolver->setPreconditioner(M_vStiffnessPreconditioner);
    }
}

void
SteadySolverOperator::setUp(const matrixEpetraPtr_Type & F,
                            const matrixEpetraPtr_Type & B,
                            const matrixEpetraPtr_Type & M_p,
                            const matrixEpetraPtr_Type & D)
{
    // we don't support this preconditioner with stabilization at the moment
    exit(1);
}


void
SteadySolverOperator::
setOptions(const Teuchos::ParameterList& solversOptions)
{
    std::shared_ptr<Teuchos::ParameterList> pressureOptions;
    pressureOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("PressureOperator")));
    setPressureMassOptions(pressureOptions);

    std::shared_ptr<Teuchos::ParameterList> momentumOptions;
    momentumOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("MomentumOperator")));
    setMomentumOptions(momentumOptions);
}

void
SteadySolverOperator::
setMomentumOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_momentumOptions = _oList;
}


void
SteadySolverOperator::
setPressureMassOptions(const parameterListPtr_Type & _oList)
{
    ASSERT_PRE(_oList.get() != 0, "oList pointer not valid");
    M_pressureOptions = _oList;
}

void
SteadySolverOperator::
updateApproximatedMomentumOperator()
{
    M_approximatedMomentumOperator->SetRowMatrix(M_F->matrixPtr());
    M_approximatedMomentumOperator->SetParameterList(*M_momentumOptions);
    if (!M_solveMomentumBelos) M_approximatedMomentumOperator->Compute();
    //
    // Epetra_Vector diag(M_F->matrixPtr()->OperatorRangeMap());
    // M_invD.reset(new Epetra_Vector(M_F->matrixPtr()->OperatorRangeMap()));
    // M_F->matrixPtr()->ExtractDiagonalCopy(diag);
    // M_invD->Reciprocal(diag);
}

void
SteadySolverOperator::
updateApproximatedPressureMassOperator()
{
    M_approximatedPressureOperator->SetRowMatrix(M_Mp->matrixPtr());
    M_approximatedPressureOperator->SetParameterList(*M_pressureOptions);
    M_approximatedPressureOperator->Compute();
}

inline
int
SteadySolverOperator::
ApplyInverse(VectorEpetra_Type const& X_velocity,
             VectorEpetra_Type const& X_pressure,
             VectorEpetra_Type & Y_velocity,
             VectorEpetra_Type & Y_pressure) const
{
    M_approximatedPressureOperator->ApplyInverse((-X_pressure).epetraVector(), Y_pressure.epetraVector());

    *M_X_velocity = X_velocity;
    *M_X_velocity -= (*M_Btranspose) * (Y_pressure);

    if (M_solveMomentumBelos)
    {
        M_linearSolver->setRightHandSide(M_X_velocity);
        M_linearSolver->solve(M_Z);
        // M_Z->epetraVector().Multiply(1.0, *M_invD, M_X_velocity->epetraVector(), 0.0);
    }
    else
        M_approximatedMomentumOperator->ApplyInverse(M_X_velocity->epetraVector(), M_Z->epetraVector());

    Y_velocity = *M_Z;

    return 0;
}

int
SteadySolverOperator::
ApplyInverse(const vector_Type& X, vector_Type& Y) const
{
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "X and Y must have the same number of vectors");

    const VectorEpetra_Type X_vectorEpetra(X, M_monolithicMap, Unique);

    // gather input values
    M_X_velocity->subset(X_vectorEpetra, M_F->map(), 0, 0);
    M_X_pressure->subset(X_vectorEpetra, M_Mp->map(), M_F->map().mapSize(), 0);

    ApplyInverse(*M_X_velocity, *M_X_pressure, *M_Y_velocity, *M_Y_pressure);

    // output vector
    VectorEpetra_Type Y_vectorEpetra(M_monolithicMap, Unique);

    // Copy the individual parts inside
    Y_vectorEpetra.subset(*M_Y_velocity, M_X_velocity->map(), 0, 0);
    Y_vectorEpetra.subset(*M_Y_pressure, M_X_pressure->map(), 0, M_X_velocity->map().mapSize());

    Y = dynamic_cast<Epetra_MultiVector&>(Y_vectorEpetra.epetraVector());

    return 0;
}

} /* end namespace Operators */

} /*end namespace */
