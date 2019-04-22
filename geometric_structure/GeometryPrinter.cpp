#include <GeometryPrinter.hpp>

namespace RedMA
{

GeometryPrinter::
GeometryPrinter()
{
}

void
GeometryPrinter::
saveToFile(TreeStructure& tree, std::string name)
{
    using namespace tinyxml2;
    typedef std::shared_ptr<TreeNode> TreeNodePtr;
    typedef std::map<std::string,std::shared_ptr<GeometricParameter> >
                                      GeometricParametersMap;

    std::queue<TreeNodePtr> nodesQueue;
    std::queue<int> outletsQueue;
    std::map<TreeNodePtr, XMLElement*> childFatherMap;
    XMLDocument outDoc;

    TreeNodePtr root = tree.getRoot();
    nodesQueue.push(root);
    outletsQueue.push(-1);

    unsigned int count = 0;
    while (nodesQueue.size() > 0)
    {
        count++;
        TreeNodePtr curNode = nodesQueue.front();
        nodesQueue.pop();
        int outletIndex = outletsQueue.front();
        outletsQueue.pop();

        if (curNode)
        {
            XMLElement* pNewElement = outDoc.NewElement("buildingblock");
            if (count == 1)
                outDoc.InsertFirstChild(pNewElement);
            else
                childFatherMap[curNode]->InsertEndChild(pNewElement);

            std::string typeName;
            if (!std::strcmp(curNode->M_block->name().c_str(),
                            "Tube"))
                typeName = "tube";
            if (!std::strcmp(curNode->M_block->name().c_str(),
                            "BifurcationSymmetric"))
                typeName = "bifurcation_symmetric";
            pNewElement->SetAttribute("type", typeName.c_str());

            if (outletIndex != -1)
                pNewElement->SetAttribute("outlet", outletIndex);

            GeometricParametersMap gpMap = curNode->M_block->getParametersMap();

            XMLElement * pElement;
            for (GeometricParametersMap::iterator it = gpMap.begin();
                 it != gpMap.end(); it++)
            {
                pElement = outDoc.NewElement(it->first.c_str());
                pElement->SetText(it->second->getValue());
                pNewElement->InsertEndChild(pElement);
            }

            // remember outlets indexing in xml starts at 1
            unsigned int outletCount = 1;
            for (auto it = curNode->M_children.begin();
                 it != curNode->M_children.end(); it++)
            {
                nodesQueue.push(*it);
                outletsQueue.push(outletCount);
                outletCount++;
                childFatherMap[*it] = pNewElement;
            }
        }

    }

    outDoc.SaveFile(name.c_str());

}


}  // namespace RedMA
