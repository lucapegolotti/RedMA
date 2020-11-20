#include "BlockAssembler.hpp"

namespace RedMA
{

template <>
void
BlockAssembler<BlockVector<VectorEp>, BlockMatrix<MatrixEp>>::
setup()
{
    typedef std::map<unsigned int, shp<TreeNode>>         NodesMap;
    typedef aAssembler<VInner, MInner>                    InnerAssembler;
    typedef std::vector<shp<TreeNode>>                    NodesVector;

    printlog(GREEN, "[BlockAssembler] initializing block assembler ... \n", this->M_data.getVerbose());

    NodesMap nodesMap = M_tree.getNodesMap();

    // allocate assemblers
    for (NodesMap::iterator it = nodesMap.begin(); it != nodesMap.end(); it++)
    {
        SHP(InnerAssembler) newAssembler;
        newAssembler = AssemblerFactory<VInner, MInner> (this->M_data, it->second);
        newAssembler->setup();
        M_primalAssemblers[it->second->M_ID] = newAssembler;
    }

    // allocate interface assemblers
    unsigned int interfaceID = 0;
    for (NodesMap::iterator it = nodesMap.begin(); it != nodesMap.end(); it++)
    {
        NodesVector children = it->second->M_children;

        unsigned int countOutlet = 0;
        unsigned int myID = it->second->M_ID;

        SHP(InnerAssembler) fatherAssembler = M_primalAssemblers[myID];

        unsigned int countChildren = 0;
        for (NodesVector::iterator itVector = children.begin();
             itVector != children.end(); itVector++)
        {
            if (*itVector)
            {
                unsigned int otherID = (*itVector)->M_ID;
                SHP(InnerAssembler) childAssembler = M_primalAssemblers[otherID];
                Interface<VInner, MInner> newInterface(fatherAssembler, myID,
                                                       childAssembler, otherID,
                                                       interfaceID);
                newInterface.M_indexOutlet = countChildren;
                SHP(InterfaceAssembler<VInner COMMA MInner>) inAssembler;
                inAssembler.reset(new InterfaceAssembler<VInner, MInner>(this->M_data,
                                                                         newInterface));
                M_dualAssemblers.push_back(inAssembler);
                interfaceID++;
            }
            countChildren++;
        }
    }

    if (!std::strcmp(this->M_data("bc_conditions/inletdirichlet", "weak").c_str(),"weak"))
    {
        SHP(InnerAssembler) inletAssembler = M_primalAssemblers[0];

        // we set the inlet to child such that we are consistent with the normal orientation
        // with respect to flow direction
        Interface<VInner, MInner> newInterface(nullptr, -1,
                                               inletAssembler, 0,
                                               interfaceID);

        SHP(InterfaceAssembler<VInner COMMA MInner>) inletInAssembler;
        inletInAssembler.reset(new InletInflowAssembler<VInner, MInner>(this->M_data,
                                                                        newInterface));

        M_dualAssemblers.push_back(inletInAssembler);
        interfaceID++;
    }

    M_numberBlocks = M_primalAssemblers.size() + M_dualAssemblers.size();

    printlog(GREEN, "done\n", this->M_data.getVerbose());
}

}
