namespace RedMA
{

template <class InVectorType, class InMatrixType>
DefaultAssemblersLibrary<InVectorType, InMatrixType>::
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
            SHP(TreeNode) defTreeNode = generateDefaultTreeNode(nameMesh);
            SHP(AssemblerType) defAssembler = AssemblerFactory<FEVECTOR COMMA FEMATRIX>(M_data, defTreeNode);
            defAssembler->initializeFEspaces();
            // defAssembler->setup();

            M_assemblersMap[mesh]= defAssembler;
            M_count++;
        }
    }
}

template <class InVectorType, class InMatrixType>
SHP(TreeNode)
DefaultAssemblersLibrary<InVectorType, InMatrixType>::
generateDefaultTreeNode(const std::string& nameMesh)
{
    if (nameMesh.find("tube") != std::string::npos)
        return generateDefaultTube(nameMesh);
    else if (nameMesh.find("bifurcation_symmetric"))
        return generateDefaultSymmetricBifurcation(nameMesh);
    else
    {
        throw new Exception("BasisGenerator: this branch must still be implemented");
    }

    return nullptr;
}

template <class InVectorType, class InMatrixType>
SHP(TreeNode)
DefaultAssemblersLibrary<InVectorType, InMatrixType>::
generateDefaultTube(const std::string& nameMesh)
{
    unsigned int diameter = std::atoi(nameMesh.substr(5,6).c_str());
    unsigned int length = std::atoi(nameMesh.substr(7,8).c_str());
    std::string refinement = nameMesh.substr(9);

    SHP(Tube) defaultTube(new Tube(M_comm, refinement, false, diameter, length, false));
    defaultTube->readMesh();

    SHP(TreeNode) treeNode(new TreeNode(defaultTube, 1234 + M_count));

    return treeNode;
}

template <class InVectorType, class InMatrixType>
SHP(TreeNode)
DefaultAssemblersLibrary<InVectorType, InMatrixType>::
generateDefaultSymmetricBifurcation(const std::string& nameMesh)
{
    unsigned int alpha = std::atoi(nameMesh.substr(13,15).c_str());
    std::string refinement = nameMesh.substr(17);

    SHP(BifurcationSymmetric) defaultBifurcation(new BifurcationSymmetric(M_comm,
                                                 refinement, false, alpha, false));
    defaultBifurcation->readMesh();

    SHP(TreeNode) treeNode(new TreeNode(defaultBifurcation, 1234 + M_count));

    return treeNode;
}

template <class InVectorType, class InMatrixType>
SHP(aAssembler<InVectorType COMMA InMatrixType>)
DefaultAssemblersLibrary<InVectorType, InMatrixType>::
getDefaultAssembler(const std::string& namemesh)
{
    return M_assemblersMap[namemesh];
}

}
