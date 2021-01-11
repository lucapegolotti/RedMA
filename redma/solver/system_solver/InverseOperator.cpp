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
setOperator(shp<LinearOperator> oper)
{
    M_invOper->setOperator(oper);
}

void
InverseOperator::
setBlockMaps(shp<BlockMaps> maps)
{
    M_maps = maps;
}

void
InverseOperator::
setPreconditioner(shp<PreconditionerOperator> prec)
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
    shp<Teuchos::ParameterList> monolithicOptions;
    monolithicOptions.reset(
        new Teuchos::ParameterList(M_solversOptions->sublist("MonolithicOperator")));
    M_pListLinSolver = monolithicOptions;
}

int
InverseOperator::
invert(const shp<aVector>& rhs, shp<aVector>& sol)
{
    shp<VECTOREPETRA> rhsEpetra = getEpetraVector(rhs, *M_maps);
    shp<VECTOREPETRA> solEpetra = getEpetraVector(sol, *M_maps);
    M_invOper->ApplyInverse(rhsEpetra->epetraVector(), solEpetra->epetraVector());

    sol = getBlockVector(solEpetra, *M_maps);

    return M_invOper->NumIter();
}


}
