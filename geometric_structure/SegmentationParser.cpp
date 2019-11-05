#include "SegmentationParser.hpp"

namespace RedMA
{

SegmentationParser::
SegmentationParser(commPtr_Type comm, std::string pthName,
                   std::string ctgrName, std::string interpolationMethod,
                   bool verbose) :
  M_comm(comm),
  M_verbose(verbose)
{
    traversePath(pthName);
    traverseSegmentation(ctgrName);
    if (!std::strcmp(interpolationMethod.c_str(),"linear"))
        linearInterpolation();
    else
        throw new Exception("Type of interpolation not implemented!");
}

SegmentationParser::
~SegmentationParser()
{

}

void
SegmentationParser::
linearInterpolation()
{
    unsigned int idContour = 0;

    unsigned int pathsize = M_path.size();
    unsigned int nContours = M_contours.size();

    M_contoursComplete.resize(pathsize);
    M_cumulativeDistance.resize(pathsize);

    Contour curContour = M_contours[idContour];
    Contour nextContour = M_contours[idContour+1];
    for (unsigned int id = 0; id < pathsize; id++)
    {
        if (id >= M_indexBegin && id < M_indexEnd)
        {
            unsigned int idCopy = id;
            double arclength = 0;
            std::vector<double> cumulativeDistance;
            cumulativeDistance.push_back(0.0);
            // first we compute the approximated arclenght from contour 1 to contour 2
            while (idCopy != nextContour.M_flag)
            {
                arclength += (M_path[idCopy+1]-M_path[idCopy]).norm();
                cumulativeDistance.push_back(arclength);
                idCopy++;
            }
            // then we linearly interpolate between the two consecutive contours
            idCopy = id;

            Vector3D axis;
            double angle;
            Matrix3D R;
            BuildingBlock::computeRotationAxisAndAngle(curContour.M_normal,
                                                       nextContour.M_normal,
                                                       axis, angle);
            do
            {
                Contour newContour;

                // relative arclength
                double s = cumulativeDistance[id-idCopy] / arclength;
                Vector3D newCenter = curContour.M_center * (1.-s) +
                                     nextContour.M_center * s;

                R = BuildingBlock::computeRotationMatrix(axis, angle * s);

                // Vector3D newNormal = curContour.M_normal * (1.-s) +
                //                      nextContour.M_normal * s;
                Vector3D newNormal = R * curContour.M_normal;
                double newRadius = curContour.M_radius * (1.-s) +
                                   nextContour.M_radius * s;

                newContour.M_center = newCenter;
                newContour.M_normal = newNormal / newNormal.norm();
                newContour.M_radius = newRadius;
                if (id == idCopy)
                    newContour.M_flag = id;
                else if (id == nextContour.M_flag)
                    newContour.M_flag = id;

                M_contoursComplete[id] = newContour;
                id++;
            } while (id != nextContour.M_flag);

            idContour++;
            curContour = M_contours[idContour];

            if (idContour+1 < nContours)
                nextContour = M_contours[idContour+1];
            // we decrement by 1 because we want to start from the former nextContour
            // at the next iteration
            id--;
        }
    }


    for (unsigned int i = M_indexBegin+1; i < M_indexEnd; i++)
    {
        M_cumulativeDistance[i] = M_cumulativeDistance[i-1] +
            (M_contoursComplete[i].M_center - M_contoursComplete[i-1].M_center).norm();
    }
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

    Vector3D prevNormal;
    bool first = true;
    do
    {
        unsigned int pathIndex = std::atoi(fContour->FirstChildElement("path_point")->
                                           Attribute("id"));

        tinyxml2::XMLElement* cPoint = fContour->FirstChildElement("contour_points")->
                                                 FirstChildElement("point");

        if (first)
        {
            M_indexBegin = pathIndex;
            first = false;
        }

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
        bool firstInner = true;
        do
        {
            Vector3D curVec = get3DVectorFromXMLElement(cPoint)-center;
            radius = radius + curVec.norm();

            if (!firstInner)
            {
                Vector3D elNormal = curVec.cross(prevVec);
                elNormal = elNormal / elNormal.norm();
                normal = normal + elNormal;
            }

            prevVec = curVec;
            firstInner = false;
            cPoint = cPoint->NextSiblingElement("point");
        } while (cPoint);

        normal = normal / normal.norm();
        if (!first)
            normal = normal.dot(prevNormal) > 0? normal : (-1.0) * normal;
        radius = radius / count;

        Contour newContour(center, normal, radius, pathIndex);

        M_contours.push_back(newContour);

        fContour = fContour->NextSiblingElement("contour");

        prevNormal = normal;
        M_indexEnd = pathIndex;
    } while (fContour);

    printlog(MAGENTA, "done\n", M_verbose);
}

TreeStructure
SegmentationParser::
createTree(const int& lengthTubes,
           const double& constVector, const double& constNormal,
           int indexBegin, int indexEnd)
{
    if (indexBegin == -1)
    {
        indexBegin = M_indexBegin;
    }

    if (indexEnd == -1)
    {
        indexEnd = M_indexEnd;
    }

    TreeStructure retTree(M_verbose);

    Vector3D prevCenter = M_contoursComplete[indexBegin].M_center;
    Vector3D prevNormal = M_contoursComplete[indexBegin].M_normal;
    double prevRadius = M_contoursComplete[indexBegin].M_radius;


    unsigned int prevIndex = M_indexBegin;
    unsigned int ind = M_indexBegin;
    unsigned int count = -1;
    // add other blocks
    double length = 0;
    int curLength = lengthTubes;
    while (1)
    {
        int status = 0;

        std::shared_ptr<Tube> childTube(new Tube(M_comm,"fine",false,1,curLength));

        Vector3D cb = prevCenter;

        status = status || childTube->setParameterValue("bx", cb[0]);
        status = status || childTube->setParameterValue("by", cb[1]);
        status = status || childTube->setParameterValue("bz", cb[2]);

        double scale = prevRadius / childTube->getInletRadius();
        status = status || childTube->setParameterValue("scale", scale);

        Vector3D axis;
        double alpha;

        BuildingBlock::computeRotationAxisAndAngle(childTube->getInletNormal(),
                                                   prevNormal, axis, alpha);

        status = status || childTube->setParameterValue("rotation_axis_x", axis[0]);
        status = status || childTube->setParameterValue("rotation_axis_y", axis[1]);
        status = status || childTube->setParameterValue("rotation_axis_z", axis[2]);
        status = status || childTube->setParameterValue("alpha", alpha);

        length += childTube->getDefLength() * scale;

        // find new target point
        while (M_cumulativeDistance[ind] - M_cumulativeDistance[indexBegin] < length
               && ind < M_indexEnd)
            ind++;

        if (ind == M_indexEnd)
        {
            if (curLength == 1)
                break;
            else
                status = 1;
        }

        // choose the closest point between ind-1 and ind
        ind = std::abs(M_cumulativeDistance[ind]-length) <
              std::abs(M_cumulativeDistance[ind-1]-length) ? ind : ind - 1;

        Vector3D targetCentral = M_contoursComplete[ind].M_center;
        Vector3D targetNormal = M_contoursComplete[ind].M_normal;
        double targetRadius = M_contoursComplete[ind].M_radius;

        double defRadius = childTube->getInletRadius();
        status = status || childTube->setParameterValue("Rout_ratio",
                                                targetRadius/(defRadius*scale));

        Matrix3D R1 = BuildingBlock::computeRotationMatrix(axis, alpha);
        Vector3D e1;
        e1[0] = 1.0; e1[1] = 0.0; e1[2] = 0.0;
        e1 = R1 * e1;

        // project target point onto plane surface
        Vector3D n = prevNormal;
        Vector3D targetCentralAbs = targetCentral - cb;
        Vector3D targetCentralProj = targetCentralAbs - targetCentralAbs.dot(n) * n;

        double alphaAxis = M_PI-std::acos(targetCentralProj.dot(e1) /
                           (e1.norm()*targetCentralProj.norm()));

        double bendAngle = 0.1;

        double f;
        if (!status)
            f = optimizeBending(alphaAxis, bendAngle, 1e7, 1e-15, R1, cb,
                                scale, childTube->getDefLength(), n,
                                M_contoursComplete[ind], constVector, constNormal);

        double err;
        err = aPosterioriCheck(alphaAxis, bendAngle, R1, cb,
                               scale, childTube->getDefLength(), n,
                               M_contoursComplete[ind]);

        std::cout << "error = " << err << std::endl;
        if (err > THRESHOLD)
            status = 1;

        status = status || childTube->setParameterValue("alpha_axis", alphaAxis);
        status = status || childTube->setParameterValue("bend", bendAngle);

        if (status == 0)
        {
            if (count == -1)
                retTree.setRoot(childTube);
            else
                retTree.addChild(count, childTube);
            prevIndex = ind;
            count++;
            curLength = lengthTubes;

            bend(alphaAxis, bendAngle, prevCenter, prevNormal, R1, cb,
                 scale, childTube->getDefLength(), n);

            prevRadius = targetRadius;
        }
        else
        {
            ind = prevIndex;
            length -= childTube->getDefLength() * scale;
            curLength--;
            std::string msg;
            if (curLength == 0)
            {
                msg = "Status was not fixed by decreasing tube length";
                throw new Exception(msg);
            }
            msg = "[SegmentationParser] Status != 0 : retrying\n";
            printlog(RED, msg, M_verbose);
        }
    }

    return retTree;
}

double
SegmentationParser::
Fbending(const double& alpha, const double& theta,
                const Matrix3D& A, const Vector3D& b,
                const double& scale,
                const double& L, const Vector3D& axis,
                const Contour& target,
                const double& const1,
                const double& const2)
{
    Vector3D newCenter;
    Vector3D newNormal;

    bend(alpha, theta, newCenter, newNormal, A, b, scale, L, axis);

    return const1 * (newCenter-target.M_center).norm() / target.M_center.norm() +
           const2 * std::abs(newNormal.dot(target.M_normal)-1.0);
}

double
SegmentationParser::
aPosterioriCheck(const double& alpha, const double& theta,
                 const Matrix3D& A, const Vector3D& b,
                 const double& scale,
                 const double& L, const Vector3D& axis,
                 const Contour& target)
{
    Vector3D newCenter;
    Vector3D newNormal;

    bend(alpha, theta, newCenter, newNormal, A, b, scale, L, axis);

    return (newCenter-target.M_center).norm() / target.M_center.norm();
}

void
SegmentationParser::
bend(const double& alpha, const double& theta,
     Vector3D& center, Vector3D& normal,
     const Matrix3D& A, const Vector3D& b,
     const double& scale,
     const double& L, const Vector3D& axis)
{

    Matrix3D rot = BuildingBlock::computeRotationMatrix(axis, alpha);

    double x = 0, y = 0, z = L;
    if (theta > 1e-5)
        Tube::bendFunctionAnalytic(x, y, z, theta, L);

    center[0] = x; center[1] = y; center[2] = z;
    center = rot * (scale * A * center) + b;

    x = -std::sin(theta); y = 0, z = std::cos(theta);
    // if (theta > 1e-5)
    //     Tube::bendFunctionAnalytic(x, y, z, theta, L);

    normal[0] = x; normal[1] = y; normal[2] = z;
    normal = -(1) * rot * A * normal;
}

double
SegmentationParser::
optimizeBending(double& alpha, double& theta,
                const double& maxIt, const double& tol,
                const Matrix3D& A, const Vector3D& b,
                const double& scale,
                const double& L, const Vector3D& axis,
                const Contour& target,
                const double& const1,
                const double& const2)
{
    printlog(MAGENTA, "[SegmentationParser] optimizing bending ...\n", M_verbose);


    double f = Fbending(alpha, theta, A, b, scale, L, axis, target, const1, const2);
    std::ostringstream streamOb1;
    streamOb1 << f;
    std::string msg = "initial loss function = ";
    msg += streamOb1.str();
    printlog(GREEN, msg + "\n", M_verbose);
    streamOb1.str("");
    streamOb1.clear();

    Vector3D newCenter, newNormal;
    bend(alpha, theta, newCenter, newNormal, A, b, scale, L, axis);
    msg = "xcenter = " + std::to_string(newCenter[0]) + "/" + std::to_string(target.M_center[0]) +
         " ycenter = " + std::to_string(newCenter[1]) + "/" + std::to_string(target.M_center[1]) +
         " zcenter = " + std::to_string(newCenter[2]) + "/" + std::to_string(target.M_center[2]);
    printlog(GREEN, msg + "\n", M_verbose);
    msg = "xnormal = " + std::to_string(newNormal[0]) + "/" + std::to_string(target.M_normal[0]) +
         " ynormal = " + std::to_string(newNormal[1]) + "/" + std::to_string(target.M_normal[1]) +
         " znormal = " + std::to_string(newNormal[2]) + "/" + std::to_string(target.M_normal[2]);
    printlog(GREEN, msg + "\n", M_verbose);


    unsigned int k = 1;
    const double eps = 1e-8;
    double incr = tol + 1;
    double oldf = 1e8;
    const double lambda = 1e-5;
    while (incr > tol && k <= maxIt)
    {
        // approximate the gradient via finite difference
        double derAlpha = (Fbending(alpha + eps, theta, A, b, scale, L,
                                   axis, target, const1, const2) - f)/eps;

        double derTheta = (Fbending(alpha, theta + eps, A, b, scale, L,
                                   axis, target, const1, const2) - f)/eps;


        alpha = alpha - lambda * derAlpha;
        theta = theta - lambda * derTheta;

        if (theta < 1e-2)
        {
            alpha += M_PI;
            theta = 0.1;
        }

        f = Fbending(alpha, theta, A, b, scale, L, axis, target, const1, const2);

        incr = std::abs((f - oldf)/f);

        if (k % 100000 == 0)
        {
            msg = "it = " + std::to_string(k);
            std::ostringstream streamOb2;
            streamOb2 << f;
            msg += " lf = " + streamOb2.str();
            streamOb2.str("");
            streamOb2.clear();
            streamOb2 << incr;
            msg += " increment = " + streamOb2.str();
            printlog(YELLOW, msg + "\n", M_verbose);
        }
        oldf = f;
        k++;
    }
    msg = "final loss function = ";
    streamOb1 << f;
    msg += streamOb1.str();
    printlog(GREEN, msg + "\n", M_verbose);


    bend(alpha, theta, newCenter, newNormal, A, b, scale, L, axis);
    msg = "xcenter = " + std::to_string(newCenter[0]) + "/" + std::to_string(target.M_center[0]) +
         " ycenter = " + std::to_string(newCenter[1]) + "/" + std::to_string(target.M_center[1]) +
         " zcenter = " + std::to_string(newCenter[2]) + "/" + std::to_string(target.M_center[2]);
    printlog(GREEN, msg + "\n", M_verbose);
    msg = "xnormal = " + std::to_string(newNormal[0]) + "/" + std::to_string(target.M_normal[0]) +
         " ynormal = " + std::to_string(newNormal[1]) + "/" + std::to_string(target.M_normal[1]) +
         " znormal = " + std::to_string(newNormal[2]) + "/" + std::to_string(target.M_normal[2]);
    printlog(GREEN, msg + "\n", M_verbose);
    return f;
}


}  // namespace RedMA
