#include <GeometryParser.hpp>

namespace ReMA
{

GeometryParser::GeometryParser(std::string fileName, commPtr_Type comm, bool verbose) :
  M_comm(comm),
  M_verbose(verbose)
{
    printlog(GREEN, "[GeometryParser] parsing " + fileName + " structure file ...\n", M_verbose);
    tinyxml2::XMLDocument doc;
    doc.LoadFile(fileName.c_str());

    std::shared_ptr<BuildingBlock> root = parseElement(doc.FirstChildElement());

    M_tree.setRoot(root);

    printlog(GREEN, "done\n", verbose);
}

std::shared_ptr<BuildingBlock> GeometryParser::parseElement(const tinyxml2::XMLElement *element)
{
    std::shared_ptr<BuildingBlock> returnBlock;

    if (!std::strcmp(element->Value(), "tube"))
    {
        returnBlock.reset(new Tube(M_comm, M_verbose));
    }

    std::map<std::string,double>& parametersMap = returnBlock->getParametersMap();

    for (std::map<std::string,double>::iterator it = parametersMap.begin();
         it != parametersMap.end(); it++)
    {
        std::string paramName = it->first;
        const tinyxml2::XMLElement* value = element->FirstChildElement(paramName.c_str());
        if (value)
        {
            parametersMap[paramName] = std::stod(value->GetText());
        }
        else
        {
            std::string warningMsg = "[GeometryParser] attribute " + paramName +
            " is not set in the datafile and will be set to default value!\n";

            printlog(YELLOW, warningMsg, M_verbose);
        }
    }

    return returnBlock;
}

}  // namespace ReMA
