#include "LinearSystemSolver.hpp"

namespace RedMA
{

template<>
BlockVector<BlockVector<VectorEp>>
LinearSystemSolver<BlockVector<VectorEp>, BlockMatrix<MatrixEp>>::
solve(BlockMatrix<BlockMatrix<MatrixEp>> matrix, BlockVector<BlockVector<VectorEp>> rhs)
{
    BlockVector<BlockVector<VectorEp>> sol;

    M_oper.reset(new LinearOperatorEp());
    M_oper->setup(matrix);

    return sol;
}

template<>
void
LinearSystemSolver<BlockVector<VectorEp>, BlockMatrix<MatrixEp>>::
setSolversOptions()
{
    std::string optionsPrec = M_data("fluid/options_preconditioner",
                                     "solversOptionsFast");
    optionsPrec += ".xml";
    M_solversOptions = Teuchos::getParametersFromXmlFile(optionsPrec);
    std::shared_ptr<Teuchos::ParameterList> monolithicOptions;
    monolithicOptions.reset(
        new Teuchos::ParameterList(M_solversOptions->sublist("MonolithicOperator")));
    M_pListLinSolver = monolithicOptions;
}


}
