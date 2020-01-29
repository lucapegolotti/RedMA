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
    M_innerPrecType = M_data("preconditioner/inner", "SIMPLE");
    if (!std::strcmp(M_innerPrecType.c_str(), "exactsolve"))
    {
    }
    else if (!std::strcmp(M_innerPrecType.c_str(), "SIMPLE"))
    {
        allocateInnerPreconditioners(A);
    }
    else
        throw new Exception("Requested inner preconditioner not implemented");

    computeAm1BT(A, BT);
}

void
SaddlePointPreconditionerEp::
allocateInnerPreconditioners(const BM& primalMatrix)
{
    unsigned int nrows = primalMatrix.nRows();

    for (unsigned int i = 0; i < nrows; i++)
    {
        SHP(NSPrec) newPrec;
        newPrec.reset(LifeV::Operators::NSPreconditionerFactory::
                      instance().createObject(M_innerPrecType));
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
computeAm1BT(const BlockMatrix<BlockMatrix<MatrixEp>>& A,
             const BlockMatrix<BlockMatrix<MatrixEp>>& BT)
{
    M_Am1BT.resize(A.nRows(), BT.nCols());

    for (unsigned int i = 0; i < A.nRows(); i++)
    {
        for (unsigned int j = 0; j < BT.nCols(); j++)
        {
            M_Am1BT.block(i,j).softCopy(computeSingleAm1BT(A.block(i,i), BT.block(i,j),i));
        }
    }


}

BlockMatrix<MatrixEp>
SaddlePointPreconditionerEp::
computeSingleAm1BT(const BlockMatrix<MatrixEp>& A, const BlockMatrix<MatrixEp>& BT,
                   const unsigned int& index)
{
    BlockMatrix<MatrixEp> retMat;
    retMat.resize(A.nRows(), BT.nCols());

    if (!std::strcmp(M_innerPrecType.c_str(),"SIMPLE"))
    {
        MAPEPETRA rangeMapU = BT.block(0,0).data()->rangeMap();
        MAPEPETRA rangeMapP = A.block(0,1).data()->rangeMap();
        MAPEPETRA domainMap = BT.block(0,0).data()->domainMap();

        VECTOREPETRA selector(domainMap);
        VECTOREPETRA colU(rangeMapU);
        VECTOREPETRA colP(rangeMapP);

        unsigned int ncols = domainMap.mapSize();

        std::vector<VectorEp> ressU(ncols);
        std::vector<VectorEp> ressP(ncols);
        for (unsigned int i = 0; i < ncols; i++)
        {
            selector.zero();
            colU.zero();
            colP.zero();
            if (domainMap.isOwned(i))
                selector[i] = 1.0;

            colU = (*BT.block(0,0).data()) * selector;

            ressU[i].data().reset(new VECTOREPETRA(rangeMapU));
            ressP[i].data().reset(new VECTOREPETRA(rangeMapP));
            M_innerPreconditioners[index]->ApplyInverse(colU, colP,
                                                        *ressU[i].data(),
                                                        *ressP[i].data());

        }
        MatrixEp Am1BTu(ressU);
        MatrixEp Am1BTp(ressP);
        retMat.block(0,0).softCopy(Am1BTu);
        retMat.block(1,0).softCopy(Am1BTp);
    }

    return retMat;
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
