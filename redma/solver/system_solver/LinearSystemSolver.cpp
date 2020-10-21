#include "LinearSystemSolver.hpp"

namespace RedMA
{

LinearSystemSolver::
LinearSystemSolver(const DataContainer& data) :
  M_data(data),
  M_numSolves(0)
{
}

void
LinearSystemSolver::
solve(const BM& matrix, const BV& rhs, BV& sol)
{
    BM matrixSparse = matrix;
    BV rhsSparse = rhs;
    if (!std::static_pointer_cast<BlockMatrix>(matrixSparse)->globalTypeIs(SPARSE))
        matrixSparse = std::static_pointer_cast<BlockMatrix>(matrixSparse)->convertInnerTo(SPARSE,M_comm);

    // std::cout << "here" << std::endl << std::flush;
    // if (!std::static_pointer_cast<BlockVector>(rhsSparse)->globalTypeIs(DISTRIBUTED))
    //     rhsSparse = std::static_pointer_cast<BlockVector>(rhsSparse)->convertInnerTo(DISTRIBUTED,std::static_pointer_cast<BlockMatrix>(matrixSparse)->commPtr());
    // std::cout << "there" << std::endl << std::flush;

    if (!M_maps)
        M_maps.reset(new BlockMaps(std::static_pointer_cast<BlockMatrix>(matrixSparse)));
    else
        M_maps->updateCollapsedMatrix(std::static_pointer_cast<BlockMatrix>(matrixSparse));

    M_oper.reset(new LinearOperator(matrixSparse,M_maps));

    M_invOper.reset(new InverseOperator(M_data));
    M_invOper->setOperator(M_oper);
    M_invOper->setBlockMaps(M_maps);

    buildPreconditioner(matrixSparse);

    M_invOper->setPreconditioner(M_prec);

    Chrono chrono;
    chrono.start();
    printlog(MAGENTA, "[LinearSystemSolver] solve ...", M_data.getVerbose());

    M_statistics.M_numIterations = M_invOper->invert(rhs, sol);

    M_statistics.M_solveTime = chrono.diff();
    std::string msg = "done, in ";
    msg += std::to_string(M_statistics.M_solveTime);
    msg += " seconds\n";
    printlog(GREEN, msg, M_data.getVerbose());

    M_numSolves++;
}

void
LinearSystemSolver::
buildPreconditioner(const BM& matrix)
{
    std::string precType = M_data("preconditioner/outer", "saddlepoint");

    if (!std::strcmp(precType.c_str(), "saddlepoint"))
    {
        unsigned int recomputeevery = M_data("preconditioner/recomputeevery", 1);
        if (M_prec == nullptr || (M_numSolves % recomputeevery) == 0)
            M_prec.reset(new SaddlePointPreconditioner(M_data, std::static_pointer_cast<BlockMatrix>(matrix)));
        else
        {
            std::static_pointer_cast<SaddlePointPreconditioner>(M_prec)->setup(std::static_pointer_cast<BlockMatrix>(matrix), false);
        }
    }
    else
    {
        throw new Exception("Unkown type of preconditioner!");
    }
    M_statistics.M_precSetupTime = M_prec->getSetupTime();
}

}
