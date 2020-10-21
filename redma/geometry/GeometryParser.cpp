#include "GeometryParser.hpp"

namespace RedMA
{

GeometryParser::
GeometryParser(const GetPot& datafile, std::string fileName,
               commPtr_Type comm, bool verbose) :
  M_comm(comm),
  M_datafile(datafile),
  M_verbose(verbose),
  M_numBlocks()
{
    printlog(MAGENTA, "[GeometryParser] parsing " +
                    fileName + " structure file ...\n", M_verbose);
    tinyxml2::XMLDocument doc;
    int status = doc.LoadFile(fileName.c_str());

    M_maxNumBlocks = M_datafile("geometric_structure/maxnumblocks", -1);

    if (status)
    {
        std::string errorMsg = "[GeometryParser] " + fileName + " does not " +
        " exist, or it is badly formatted!";
        throw Exception(errorMsg);
    }

    tinyxml2::XMLElement* rootElement = doc.FirstChildElement();
    traverseXML(rootElement, 0);

    printlog(MAGENTA, "done\n", verbose);
}

void
GeometryParser::
traverseXML(tinyxml2::XMLElement* curElement,
            unsigned int IDfather)
{
    if (M_maxNumBlocks < 0 || M_numBlocks < M_maxNumBlocks)
    {
        std::shared_ptr<BuildingBlock> newNode;
        unsigned int outletParent = -1;
        if (curElement)
            newNode = parseElement(curElement, outletParent);
        else
            return;

        unsigned int newID = 0;
        if (M_tree.isEmpty())
            M_tree.setRoot(newNode);
        else
            newID = M_tree.addChild(IDfather, newNode, outletParent);

        if (curElement->NoChildren())
            return;

        tinyxml2::XMLElement* childElement =
                              curElement->FirstChildElement();
        while(childElement)
        {
            if (childElement->Attribute("type"))
                traverseXML(childElement, newID);
            childElement = childElement->NextSiblingElement();
        }
    }
}

GeometryParser::BuildingBlockPtr
GeometryParser::
parseElement(const XMLEl *element, unsigned int& outletParent)
{
    BuildingBlockPtr returnBlock;

    outletParent = -1;
    if (element->Attribute("outlet"))
    {
        // decrease by one because we want numbering of outlets to start from
        // 1
        outletParent = std::stoi(element->Attribute("outlet")) - 1;
    }

    std::string ref = "coarse";
    if (element->Attribute("refinement"))
    {
        ref = element->Attribute("refinement");
    }

    if (!std::strcmp(element->Attribute("type"), "tube"))
    {
        printlog(GREEN, std::string("[GeometryParser] parsing building block") +
                       " of type tube\n", M_verbose);

        // diameter
        int d = 1;
        if (element->Attribute("d"))
        {
            d = std::atoi(element->Attribute("d"));
        }

        // length
        int L = 1;
        if (element->Attribute("L"))
        {
            L = std::atoi(element->Attribute("L"));
        }

        returnBlock.reset(new Tube(M_comm, ref, M_verbose, d, L));
    }
    else if (!std::strcmp(element->Attribute("type"),
                          "bifurcation_symmetric"))
    {
        std::string msg = std::string("[GeometryParser] parsing") +
                  " building block of type bifurcation symmetric\n";
        printlog(GREEN, msg, M_verbose);

        // angle
        int angle = 50;
        if (element->Attribute("angle"))
            angle = std::atoi(element->Attribute("angle"));

        returnBlock.reset(new BifurcationSymmetric(M_comm, ref, M_verbose, angle));
    }
    else if (!std::strcmp(element->Attribute("type"),
                          "aorta"))
    {
        std::string msg = std::string("[GeometryParser] parsing") +
                  " building block of type aorta\n";
        printlog(GREEN, msg, M_verbose);

        returnBlock.reset(new Aorta(M_comm, "aorta", M_verbose));
    }
    else if (!std::strcmp(element->Attribute("type"),
                          "aortabif0"))
    {
        std::string msg = std::string("[GeometryParser] parsing") +
                  " building block of type aortabif0\n";
        printlog(GREEN, msg, M_verbose);

        returnBlock.reset(new AortaBifurcation0(M_comm, ref, "aortabif0", M_verbose));
    }
    else if (!std::strcmp(element->Attribute("type"),
                          "aortabif1"))
    {
        std::string msg = std::string("[GeometryParser] parsing") +
                  " building block of type aortabif0\n";
        printlog(GREEN, msg, M_verbose);

        returnBlock.reset(new AortaBifurcation1(M_comm, ref, "aortabif1", M_verbose));
    }
    else
    {
        std::string warningMsg = "[GeometryParser] building block "
        + std::string(element->Attribute("type")) +
        " is not implemented! Skipping\n";

        printlog(YELLOW, warningMsg, M_verbose);
    }

    std::string discrMethod = "fem";
    if (element->Attribute("method"))
    {
        discrMethod = element->Attribute("method");
    }

    std::string assembler = "stokes";
    if (element->Attribute("assembler"))
    {
        assembler = element->Attribute("assembler");
    }

    returnBlock->setDiscretizationMethod(discrMethod);
    returnBlock->setAssemblerType(assembler);

    typedef std::shared_ptr<GeometricParameter> GeometricParameterPtr;
    std::map<std::string,GeometricParameterPtr>& parametersMap =
                                  returnBlock->getParametersMap();

    typedef std::map<std::string,GeometricParameterPtr>::iterator MapIterator;

    for (MapIterator it = parametersMap.begin();
         it != parametersMap.end(); it++)
    {
        std::string paramName = it->first;
        const tinyxml2::XMLElement* value =
                    element->FirstChildElement(paramName.c_str());

        if (value)
            *parametersMap[paramName] = std::stod(value->GetText());
    }
    returnBlock->setDatafile(M_datafile);

    M_numBlocks++;

    return returnBlock;
}

TreeStructure&
GeometryParser::
getTree()
{
    return M_tree;
}

}  // namespace RedMA
