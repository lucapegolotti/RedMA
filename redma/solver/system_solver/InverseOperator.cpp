#include "InverseOperator.hpp"

namespace RedMA
{

InverseOperator::
InverseOperator(const DataContainer& data) :
  M_data(data)
{
    using namespace LifeV::Operators;
    setSolverOptions();

    std::string solverType(M_pListLinSolver->get<std::string>("Linear Solver Type"));
    M_invOper.reset(InvertibleOperatorFactory::instance().createObject(solverType));
    M_invOper->setParameterList(M_pListLinSolver->sublist(solverType));
}


void
InverseOperator::
setOperator(SHP(LinearOperator) oper)
{
    // M_maps = oper->getBlockMaps();
    M_invOper->setOperator(oper);
}

void
InverseOperator::
setPreconditioner(SHP(PreconditionerOperator) prec)
{
    M_invOper->setPreconditioner(prec);
}

void
InverseOperator::
setSolverOptions()
{
    std::string optionsPrec = M_data("preconditioner/options",
                                     "datafiles/solversOptionsFast");
    optionsPrec += ".xml";
    M_solversOptions = Teuchos::getParametersFromXmlFile(optionsPrec);
    std::shared_ptr<Teuchos::ParameterList> monolithicOptions;
    monolithicOptions.reset(
        new Teuchos::ParameterList(M_solversOptions->sublist("MonolithicOperator")));
    M_pListLinSolver = monolithicOptions;
}

int
InverseOperator::
invert(const SHP(aVector)& rhs, SHP(aVector)& sol)
{
    SHP(VECTOREPETRA) rhsEpetra = getEpetraVector(rhs, *M_maps);
    SHP(VECTOREPETRA) solEpetra = getEpetraVector(sol, *M_maps);

    M_invOper->ApplyInverse(rhsEpetra->epetraVector(), solEpetra->epetraVector());

    sol = getBlockVector(solEpetra, *M_maps);

    return M_invOper->NumIter();
}


}
