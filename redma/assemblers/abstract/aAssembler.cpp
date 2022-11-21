#include "aAssembler.hpp"

namespace RedMA
{

aAssembler::
aAssembler(const DataContainer& data) :
  M_data(data),
  M_nComponents(0),
  M_name("aAssembler")
{
}

aAssembler::
aAssembler(const DataContainer& data,
           shp<TreeNode> treeNode) :
  M_data(data),
  M_treeNode(treeNode),
  M_comm(treeNode->M_block->getComm()),
  M_nComponents(0),
  M_name("aAssembler")
{
}

}
