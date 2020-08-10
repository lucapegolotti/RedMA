#include "AssemblerFactory.hpp"

namespace RedMA
{

SHP(aAssembler)
AssemblerFactory(const DataContainer& data, SHP(TreeNode) treeNode)
{
    std::shared_ptr<aAssembler> ret;
    std::string assemblerString = data("assembler/type","stokes");

    // if (!std::strcmp(assemblerString.c_str(),"stokes"))
    //     ret.reset(new StokesAssembler<InVectorType,InMatrixType>(data,
    //                                                              treeNode));
    // else if (!std::strcmp(assemblerString.c_str(),"navierstokes"))
    //     ret.reset(new NavierStokesAssembler<InVectorType,InMatrixType>(data,
    //                                                                    treeNode));
    // else
    //     throw new Exception("Specified assembler is not implemented!");

    return ret;
}

}
