#include "DefaultAssemblersLibrary.hpp"
#include <redma/assemblers/abstract/aAssembler.hpp>

namespace RedMA
{

DefaultAssemblersLibrary::
DefaultAssemblersLibrary(const DataContainer& data, const std::set<std::string>& meshes, EPETRACOMM comm) :
  M_data(data),
  M_comm(comm),
  M_count(0)
{
    for (auto mesh : meshes)
    {
        unsigned int dashpos = mesh.find_last_of("/");
        unsigned int formatpos = mesh.find_last_of(".");
        std::string nameMesh = mesh.substr(dashpos + 1,mesh.size()+dashpos+1-formatpos);

        if (M_assemblersMap.find(mesh) == M_assemblersMap.end())
        {
            shp<TreeNode> defTreeNode = generateDefaultTreeNode(nameMesh);
            if (defTreeNode)
            {
                // we set these to stokes and fem because we are only interested in
                // the finite element spaces
                defTreeNode->M_block->setAssemblerType("stokes");
                defTreeNode->M_block->setDiscretizationMethod("fem");
                shp<aAssembler> defAssembler = AssemblerFactory(M_data, defTreeNode);
                defAssembler->initializeFEspaces();
                M_assemblersMap[mesh]= defAssembler;
                M_count++;
            }
        }
    }
}

shp<TreeNode>
DefaultAssemblersLibrary::
generateDefaultTreeNode(const std::string& nameMesh)
{
    if (nameMesh.find("tube") != std::string::npos)
        return generateDefaultTube(nameMesh);
    else if (nameMesh.find("bif_sym") != std::string::npos)
        return generateDefaultSymmetricBifurcation(nameMesh);
    else if (nameMesh.find("bypass") != std::string::npos)
        return generateDefaultBypass(nameMesh);
    else
        printlog(YELLOW, "[DefaultAssemblersLibrary]: no default assembler available for " + nameMesh + " mesh!", true);

    return nullptr;
}

shp<TreeNode>
DefaultAssemblersLibrary::
generateDefaultTube(const std::string& nameMesh)
{
    unsigned int diameter = std::atoi(nameMesh.substr(5,6).c_str());
    unsigned int length = std::atoi(nameMesh.substr(7,8).c_str());
    std::string refinement = nameMesh.substr(9);

    shp<Tube> defaultTube(new Tube(M_comm, refinement, false, diameter, length, false));
    defaultTube->readMesh();

    shp<TreeNode> treeNode(new TreeNode(defaultTube, 1234 + M_count));

    return treeNode;
}

shp<TreeNode>
DefaultAssemblersLibrary::
generateDefaultSymmetricBifurcation(const std::string& nameMesh)
{
    unsigned int alpha = std::atoi(nameMesh.substr(13,15).c_str());
    std::string refinement = nameMesh.substr(17);

    shp<BifurcationSymmetric> defaultBifurcation(new BifurcationSymmetric(M_comm,
                                                 refinement, false, alpha, false));
    defaultBifurcation->readMesh();

    shp<TreeNode> treeNode(new TreeNode(defaultBifurcation, 1234 + M_count));

    return treeNode;
}

shp<TreeNode>
DefaultAssemblersLibrary::
generateDefaultBypass(const std::string &nameMesh)
{
    shp<Bypass> defaultBypass(new Bypass(M_comm, "coarse", false));

    defaultBypass->readMesh();

    shp<TreeNode> treeNode(new TreeNode(defaultBypass, 1234 + M_count));

    return treeNode;
}

shp<aAssembler>
DefaultAssemblersLibrary::
getDefaultAssembler(const std::string& namemesh)
{
    if (M_assemblersMap.find(namemesh) == M_assemblersMap.end())
        return nullptr;
    return M_assemblersMap[namemesh];
}

}
