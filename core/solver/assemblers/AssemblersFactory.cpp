#include <AssemblersFactory.hpp>

namespace RedMA
{

std::shared_ptr<AbstractAssembler>
AssemblersFactory(const GetPot& datafile, std::shared_ptr<Epetra_Comm> comm,
                  std::shared_ptr<TreeNode> treenode, bool verbose)
{
    std::shared_ptr<AbstractAssembler> newAssembler;
    std::string assemblertype = datafile("assembler/type","navierstokes");

    if (!std::strcmp(assemblertype.c_str(),"navierstokes"))
        newAssembler = std::shared_ptr<NavierStokesAssembler>
                    (new NavierStokesAssembler(datafile,comm,treenode,verbose));
    else if (!std::strcmp(assemblertype.c_str(),"pseusofsi"))
        newAssembler = std::shared_ptr<PseudoFSIAssembler>
                    (new PseudoFSIAssembler(datafile,comm,treenode,verbose));

    return newAssembler;
}

} // namespace RedMA
