#include "SegmentationParser.hpp"

namespace RedMA
{

SegmentationParser::
SegmentationParser(std::string pthName, std::string ctgrName, bool verbose) :
  M_verbose(verbose)
{
    traversePath(pthName);
    traverseSegmentation(ctgrName);
}

SegmentationParser::
~SegmentationParser()
{

}

SegmentationParser::Vector3D
SegmentationParser::
get3DVectorFromXMLElement(tinyxml2::XMLElement* data)
{
    std::string xs;
    if (data->Attribute("x"))
        xs = data->Attribute("x");

    std::string ys;
    if (data->Attribute("y"))
        ys = data->Attribute("y");

    std::string zs;
    if (data->Attribute("z"))
        zs = data->Attribute("z");

    double x = std::stod(xs);
    double y = std::stod(ys);
    double z = std::stod(zs);
    Vector3D newVector(x,y,z);
    return newVector;
}

void
SegmentationParser::
traversePath(std::string pthName)
{
    printlog(MAGENTA, "[SegmentationParser] parsing " +
                    pthName + " path file ...\n", M_verbose);
    tinyxml2::XMLDocument doc;
    int status = doc.LoadFile(pthName.c_str());

    if (status)
    {
        std::string errorMsg = "[SegmentationParser] " + pthName + " does not " +
        " exist, or it is badly formatted!";
        throw Exception(errorMsg);
    }

    tinyxml2::XMLElement* rootElement = doc.FirstChildElement("path");

    tinyxml2::XMLElement* pPoints = rootElement->FirstChildElement("timestep")->
                                                 FirstChildElement("path_element")->
                                                 FirstChildElement("path_points");

    tinyxml2::XMLElement* fPoint = pPoints->FirstChildElement("path_point");

    std::vector<std::string> tags;
    tags.push_back("pos");
    tags.push_back("tangent");
    tags.push_back("rotation");

    std::vector<Vector3D>* vectors[3];
    vectors[0] = &M_path;
    vectors[1] = &M_tangents;
    vectors[2] = &M_rotation;

    do
    {
        unsigned int ind = 0;
        for (auto it = tags.begin(); it != tags.end(); it++)
        {
            tinyxml2::XMLElement* data = fPoint->FirstChildElement(it->c_str());

            vectors[ind]->push_back(get3DVectorFromXMLElement(data));
            ind++;
        }

        fPoint = fPoint->NextSiblingElement("path_point");
    } while (fPoint);

    printlog(MAGENTA, "done\n", M_verbose);
}

void
SegmentationParser::
traverseSegmentation(std::string ctgrName)
{
    printlog(MAGENTA, "[SegmentationParser] parsing " +
                    ctgrName + " segmentation file ...\n", M_verbose);
    tinyxml2::XMLDocument doc;
    int status = doc.LoadFile(ctgrName.c_str());

    if (status)
    {
        std::string errorMsg = "[SegmentationParser] " + ctgrName + " does not " +
        " exist, or it is badly formatted!";
        throw Exception(errorMsg);
    }

    tinyxml2::XMLElement* rootElement = doc.FirstChildElement("contourgroup");
    tinyxml2::XMLElement* pTimestep = rootElement->FirstChildElement("timestep");
    tinyxml2::XMLElement* fContour = pTimestep->FirstChildElement("contour");
    do
    {
        unsigned int pathIndex = std::atoi(fContour->FirstChildElement("path_point")->
                                           Attribute("id"));

        tinyxml2::XMLElement* cPoint = fContour->FirstChildElement("contour_points")->
                                                 FirstChildElement("point");

        // first find center of the face
        Vector3D center;
        unsigned int count = 0;
        do
        {
            Vector3D curVec = get3DVectorFromXMLElement(cPoint);
            center = center + curVec;

            cPoint = cPoint->NextSiblingElement("point");
            count++;
        } while (cPoint);

        center = center / count;

        // now we search for the normal and radius
        cPoint = fContour->FirstChildElement("contour_points")->
                           FirstChildElement("point");

        Vector3D normal(0, 0, 0);
        Vector3D prevVec(0, 0, 0);
        double radius = 0;
        bool first = true;
        do
        {
            Vector3D curVec = get3DVectorFromXMLElement(cPoint)-center;
            radius = radius + curVec.norm();

            if (!first)
            {
                Vector3D elNormal = curVec.cross(prevVec);
                elNormal = elNormal / elNormal.norm();
                normal = normal + elNormal;
            }

            prevVec = curVec;
            first = false;
            cPoint = cPoint->NextSiblingElement("point");
        } while (cPoint);

        normal = normal / normal.norm();
        radius = radius / count;

        Contour newContour(center, normal, radius, pathIndex);

        newContour.print();

        M_contours.push_back(newContour);

        fContour = fContour->NextSiblingElement("contour");
    } while (fContour);

    printlog(MAGENTA, "done\n", M_verbose);
}

}  // namespace RedMA
