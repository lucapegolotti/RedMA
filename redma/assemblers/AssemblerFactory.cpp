#include "AssemblerFactory.hpp"

namespace RedMA
{

SHP(aAssembler)
AssemblerFactory(const DataContainer& data, SHP(TreeNode) treeNode)
{
    std::shared_ptr<aAssembler> ret;
    std::string assemblerString = data("assembler/type","stokes");

    std::string method = treeNode->M_block->getDiscretizationMethod();
    std::string assemblerType = treeNode->M_block->getAssemblerType();

    // if (!std::strcmp(assemblerString.c_str(),"stokes"))
    //     ret.reset(new StokesAssembler<InVectorType,InMatrixType>(data,
    //                                                              treeNode));
    // else if (!std::strcmp(assemblerString.c_str(),"navierstokes"))
    //     ret.reset(new NavierStokesAssembler<InVectorType,InMatrixType>(data,
    //                                                                    treeNode));

    if (!std::strcmp(assemblerType.c_str(),"stokes"))
    {
        if (!std::strcmp(method.c_str(),"fem"))
        {
            ret.reset(new StokesAssemblerFE(data, treeNode));
        }
        else if (!std::strcmp(method.c_str(),"rb"))
        {
            ret.reset(new StokesAssemblerRB(data, treeNode));
        }
    }
    else if (!std::strcmp(assemblerType.c_str(),"navierstokes"))
    {
        if (!std::strcmp(method.c_str(),"fem"))
        {
            ret.reset(new NavierStokesAssemblerFE(data, treeNode));
        }
        else if (!std::strcmp(method.c_str(),"rb"))
        {
            ret.reset(new NavierStokesAssemblerRB(data, treeNode));
        }
    }
    else if (!std::strcmp(assemblerType.c_str(),"navierstokessupg"))
    {
        if (!std::strcmp(method.c_str(),"fem"))
        {
            ret.reset(new NavierStokesAssemblerFE(data, treeNode, "supg"));
        }
    }
    else
        throw new Exception("Specified assembler is not implemented!");

    return ret;
}

}
