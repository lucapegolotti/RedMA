#include <GeometryParser.hpp>

namespace ReMA
{

GeometryParser::
GeometryParser(std::string fileName, commPtr_Type comm,
               bool verbose) :
  M_comm(comm),
  M_verbose(verbose)
{
    printlog(GREEN, "[GeometryParser] parsing " +
                    fileName + " structure file ...\n", M_verbose);
    tinyxml2::XMLDocument doc;
    doc.LoadFile(fileName.c_str());

    tinyxml2::XMLElement* rootElement = doc.FirstChildElement();

    traverseXML(rootElement,0);

    printlog(GREEN, "done\n", verbose);
}

void
GeometryParser::
traverseXML(tinyxml2::XMLElement* curElement,
            unsigned int IDfather)
{
    std::shared_ptr<BuildingBlock> newNode;
    if (curElement)
        newNode = parseElement(curElement);
    else
        return;

    unsigned int newID = 0;
    if (M_tree.isEmpty())
    {
        M_tree.setRoot(newNode);
    }
    else
    {
        newID = M_tree.addChild(IDfather, newNode);
    }

    if (curElement->NoChildren())
        return;

    tinyxml2::XMLElement* childElement =
                          curElement->FirstChildElement();
    while(childElement)
    {
        if (childElement->Attribute("type"))
        {
            traverseXML(childElement, newID);
        }
        childElement = childElement->NextSiblingElement();
    }
}

GeometryParser::BuildingBlockPtr
GeometryParser::
parseElement(const XMLEl *element)
{
    BuildingBlockPtr returnBlock;

    if (!std::strcmp(element->Attribute("type"), "tube"))
    {
        printlog(CYAN, "[GeometryParser] parsing building block " +
                       std::string("of type tube)\n", M_verbose));
        returnBlock.reset(new Tube(M_comm, M_verbose));
    }
    else if (!std::strcmp(element->Attribute("type"),
                          "bifurcation_symmetric"))
    {
        std::string msg = std::string("[GeometryParser] parsing") +
                  "building block of type bifurcation symmetric\n";
        printlog(CYAN, msg, M_verbose);
        returnBlock.reset(new BifurcationSymmetric(M_comm, M_verbose));
    }
    else
    {
        std::string warningMsg = "[GeometryParser] building block "
        + std::string(element->Attribute("type")) +
        " is not implemented! Skipping\n";

        printlog(YELLOW, warningMsg, M_verbose);
    }

    std::map<std::string,double>& parametersMap =
                                  returnBlock->getParametersMap();

    typedef std::map<std::string,double>::iterator MapIterator;

    for (MapIterator it = parametersMap.begin();
         it != parametersMap.end(); it++)
    {
        std::string paramName = it->first;
        const tinyxml2::XMLElement* value =
                    element->FirstChildElement(paramName.c_str());

        if (value)
            parametersMap[paramName] = std::stod(value->GetText());
        else
        {
            std::string warningMsg = "[GeometryParser] parameter " +
            paramName + " is not set in the datafile and will be " +
            "set to default value!\n";

            printlog(YELLOW, warningMsg, M_verbose);
        }
    }

    return returnBlock;
}

TreeStructure&
GeometryParser::
getTree()
{
    return M_tree;
}

}  // namespace ReMA
