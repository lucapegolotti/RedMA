#include <GeometryParser.hpp>

namespace ReMA
{

GeometryParser::GeometryParser(std::string fileName)
{
    tinyxml2::XMLDocument doc;
    doc.LoadFile(filename);

    

}

}  // namespace ReMA
