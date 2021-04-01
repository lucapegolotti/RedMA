#include "AssemblerFactory.hpp"

namespace RedMA
{

shp<aAssembler>
AssemblerFactory(const DataContainer& data, shp<TreeNode> treeNode)
{
    shp<aAssembler> ret;

    std::string method = treeNode->M_block->getDiscretizationMethod();
    std::string assemblerType = treeNode->M_block->getAssemblerType();

    if (!std::strcmp(assemblerType.c_str(), "stokes"))
    {
        if (!std::strcmp(method.c_str(), "fem"))
        {
            ret.reset(new StokesAssemblerFE(data, treeNode));
        }
        else if (!std::strcmp(method.c_str(), "rb"))
        {
            ret.reset(new StokesAssemblerRB(data, treeNode));
        }
    }
    else if (!std::strcmp(assemblerType.c_str(), "navierstokes"))
    {
        if (!std::strcmp(method.c_str(), "fem"))
        {
            ret.reset(new NavierStokesAssemblerFE(data, treeNode));
        }
        else if (!std::strcmp(method.c_str(), "rb"))
        {
            ret.reset(new NavierStokesAssemblerRB(data, treeNode));
        }
    }
    else if (!std::strcmp(assemblerType.c_str(), "navierstokessupg"))
    {
        if (!std::strcmp(method.c_str(), "fem"))
        {
            ret.reset(new NavierStokesAssemblerFE(data, treeNode, "supg"));
        }
        else if (!std::strcmp(method.c_str(), "rb"))
        {
            ret.reset(new NavierStokesAssemblerRB(data, treeNode, "supg"));
        }
    }
    else if (!std::strcmp(assemblerType.c_str(), "navierstokes_membrane"))
    {
        if (!std::strcmp(method.c_str(), "fem"))
        {
            ret.reset(new MembraneAssemblerFE(data, treeNode));
        }
        else if (!std::strcmp(method.c_str(), "rb"))
        {
            ret.reset(new MembraneAssemblerRB(data, treeNode));
        }
    }
    else if (!std::strcmp(assemblerType.c_str(), "navierstokessupg_membrane"))
    {
        if (!std::strcmp(method.c_str(), "fem"))
        {
            ret.reset(new MembraneAssemblerFE(data, treeNode, "supg"));
        }
        else if (!std::strcmp(method.c_str(), "rb"))
        {
            ret.reset(new MembraneAssemblerRB(data, treeNode, "supg"));
        }
    }
    else
        throw new Exception("Specified assembler is not implemented!");

    return ret;
}

}
