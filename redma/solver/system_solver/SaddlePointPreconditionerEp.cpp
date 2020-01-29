#include "SaddlePointPreconditionerEp.hpp"

namespace RedMA
{

SaddlePointPreconditionerEp::
SaddlePointPreconditionerEp(const DataContainer& data, const BM& matrix) :
  M_data(data),
  M_matrix(matrix)
{
    // split matrices
    unsigned int nBlocks = matrix.nRows();

    unsigned int nPrimal = 0;
    while (!matrix.block(nPrimal,nPrimal).isNull())
        nPrimal++;

    unsigned int nDual = nBlocks - nPrimal;

    BM A = matrix.getSubmatrix(0, nPrimal-1, 0, nPrimal-1);
    BM BT = matrix.getSubmatrix(0, nPrimal-1, nPrimal, nBlocks-1);
    BM B = matrix.getSubmatrix(nPrimal, nBlocks-1, 0, nPrimal-1);

    // read options
    setSolverOptions();

    // build Schur complement
    std::string innerPrecType = M_data("preconditioner/inner", "SIMPLE");
    if (!std::strcmp(innerPrecType.c_str(), "exactsolve"))
    {
    }
    else if (!std::strcmp(innerPrecType.c_str(), "SIMPLE"))
    {
        allocateInnerPreconditioners(A, innerPrecType);
    }
    else
    {
        throw new Exception("Requested inner preconditioner not implemented");
    }
}

void
SaddlePointPreconditionerEp::
allocateInnerPreconditioners(const BM& primalMatrix, std::string precType)
{
    unsigned int nrows = primalMatrix.nRows();

    for (unsigned int i = 0; i < nrows; i++)
    {
        SHP(NSPrec) newPrec;
        newPrec.reset(LifeV::Operators::NSPreconditionerFactory::
                      instance().createObject(precType));
        newPrec->setOptions(*M_solversOptionsInner);

        if (primalMatrix.block(i,i).block(1,1).data())
            newPrec->setUp(primalMatrix.block(i,i).block(0,0).data(),
                           primalMatrix.block(i,i).block(1,0).data(),
                           primalMatrix.block(i,i).block(0,1).data(),
                           primalMatrix.block(i,i).block(1,1).data());
        else
            newPrec->setUp(primalMatrix.block(i,i).block(0,0).data(),
                           primalMatrix.block(i,i).block(1,0).data(),
                           primalMatrix.block(i,i).block(0,1).data());

        std::vector<SHP(Epetra_Map)> localRangeMaps(2);
        std::vector<SHP(Epetra_Map)> localDomainMaps(2);

        for (unsigned int j = 0; j < 2; j++)
            localRangeMaps[j].reset(new Epetra_Map(primalMatrix.block(i,i).
                                    block(j,0).data()->matrixPtr()->OperatorRangeMap()));

        for (unsigned int j = 0; j < 2; j++)
            localDomainMaps[j].reset(new Epetra_Map(primalMatrix.block(i,i).
                                     block(0,j).data()->matrixPtr()->OperatorDomainMap()));

        SHP(LifeV::BlockEpetra_Map) rmblock(new LifeV::BlockEpetra_Map(localRangeMaps));
        SHP(LifeV::BlockEpetra_Map) dmblock(new LifeV::BlockEpetra_Map(localDomainMaps));
        newPrec->setRangeMap(rmblock);
        newPrec->setDomainMap(dmblock);

        newPrec->updateApproximatedMomentumOperator();
        newPrec->updateApproximatedSchurComplementOperator();

        M_innerPreconditioners.push_back(newPrec);
    }
}

void
SaddlePointPreconditionerEp::
setSolverOptions()
{
    std::string optionsPrec = M_data("preconditioner/options",
                                     "datafiles/solversOptionsFast");
    optionsPrec += ".xml";
    M_solversOptionsInner = Teuchos::getParametersFromXmlFile(optionsPrec);
    std::shared_ptr<Teuchos::ParameterList> monolithicOptions;
    monolithicOptions.reset(
        new Teuchos::ParameterList(M_solversOptionsInner->sublist("InnerBlockOperator")));
    M_pListLinSolver = monolithicOptions;
}


int
SaddlePointPreconditionerEp::
ApplyInverse(const super::vector_Type& X, super::vector_Type& Y) const
{
    Y = X;

    return 0;
}

}
