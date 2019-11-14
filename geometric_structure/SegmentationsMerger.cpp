#include <SegmentationsMerger.hpp>

namespace RedMA
{

SegmentationsMerger::
SegmentationsMerger(commPtr_Type comm, bool verbose) :
  M_comm(comm),
  M_verbose(verbose)
{
}


TreeStructure
SegmentationsMerger::
merge(std::vector<SegmentationParserPtr> parsers)
{
    bool** connectivity = getConnectivity();

    unsigned int npaths = parsers.size();

    TreeStructure tree;

    for (unsigned int i = 0; i < npaths; i++)
    {
        for (unsigned int j = 0; j < npaths; j++)
        {
            if (connectivity[i][j])
                mergeTwoSegmentations(parsers[i], parsers[j], tree);
        }
    }

    deallocateConnectivity(connectivity, npaths);

    return tree;
}

bool**
SegmentationsMerger::
getConnectivity()
{
    bool** connectivity = new bool*[2];

    for (unsigned int i = 0; i < 2; i++)
    {
        connectivity[i] = new bool[2];
        for (unsigned int j = 0; j < 2; j++)
        {
            connectivity[i][j] = false;
        }
    }

    connectivity[0][1] = true;

    return connectivity;
}

void
SegmentationsMerger::
deallocateConnectivity(bool** connectivity, unsigned int size)
{
    for (unsigned int i = 0; i < size; i++)
    {
        delete[] connectivity[i];
    }

    delete[] connectivity;
}

void
SegmentationsMerger::
mergeTwoSegmentations(SegmentationParserPtr segmentationFather,
                      SegmentationParserPtr segmentationChild,
                      TreeStructure& outputTree)
{
    typedef std::vector<Contour>::iterator       VectorContourIt;

    // find guess for the intersection point
    std::vector<Contour> contours1 = segmentationFather->getContours();
    std::vector<Contour> contours2 = segmentationChild->getContours();

    Vector3D contourCenter;
    double dist = 0.0;

    Vector3D firstCenter = contours2[segmentationChild->getIndexBegin()].M_center;
    double radius = contours2[segmentationChild->getIndexBegin()].M_radius;
    unsigned int startIndex = 0;
    unsigned int index = segmentationChild->getIndexBegin();
    unsigned int closestPoint;
    dist = findClosestPoint(contourCenter, segmentationFather,
                            startIndex, closestPoint);
    do
    {
        contourCenter = contours2[index].M_center;
        dist = findClosestPoint(contourCenter, segmentationFather,
                                startIndex, closestPoint);
        startIndex = closestPoint;
        index++;
    } while(dist < 0.7 * radius);

    Vector3D guessCenter = (contourCenter + contours1[closestPoint].M_center)/2;
    Vector3D guessNormal = (contours2[index-1].M_normal + contours1[closestPoint].M_normal)/2;
    // normal to the plane where the bifurcation lies
    Vector3D guessTransv = (contourCenter - firstCenter).cross(
                            contours1[closestPoint].M_center - firstCenter);
    guessTransv = guessTransv / guessTransv.norm();

    std::shared_ptr<BifurcationSymmetric> bifurcation;
    bifurcation.reset(new BifurcationSymmetric(M_comm,"coarse",false,50));

    placeBifurcation(guessCenter, guessNormal, guessTransv,
                     contours1[closestPoint].M_radius,
                     segmentationFather, segmentationChild, bifurcation);


    // forward parent until bifurcation
    Contour inletBifurcation = bifurcation->getInlet();
    TreeStructure inBranch = segmentationFather->createTreeForward(1, 1.0, 1.0,
                                                                   nullptr,
                                                            &inletBifurcation);

    unsigned int idLastElement = inBranch.getNodesMap().size()-1;

    bifurcation->resetInletOutlets();
    inBranch.addChild(idLastElement,bifurcation);

    inBranch.traverseAndDeformGeometries(false);
    Contour coutlet = inBranch.getNodesMap()[idLastElement]->M_block->getOutlet(0);
    inBranch.resetInletOutlets();

    rotateBifurcation(coutlet.M_center, coutlet.M_normal, guessTransv, coutlet.M_radius,
                      segmentationFather, segmentationChild, bifurcation);
    outputTree = outputTree + inBranch;

    Contour outlet1 = bifurcation->getOutlet(0);
    Contour outlet2 = bifurcation->getOutlet(1);
    double dist1;
    double dist2;
    dist1 = findClosestPoint(outlet1.M_center,
                             segmentationFather, 0, closestPoint);

    dist2 = findClosestPoint(outlet1.M_center,
                             segmentationChild, 0, closestPoint);

    outlet1.M_normal = (-1) * outlet1.M_normal;
    outlet2.M_normal = (-1) * outlet2.M_normal;

    TreeStructure outBranch;
    TreeStructure otherBranch;
    if (dist1 < dist2)
    {
        outBranch = segmentationFather->createTreeForward(1, 1.0, 1.0, &outlet1, nullptr);
        otherBranch = segmentationChild->createTreeForward(1, 1.0, 1.0, &outlet2, nullptr);
    }
    else
    {
        outBranch = segmentationChild->createTreeForward(1, 1.0, 1.0, &outlet1, nullptr);
        otherBranch = segmentationFather->createTreeForward(1, 1.0, 1.0, &outlet2, nullptr);
    }

    outBranch.traverseAndDeformGeometries(false);
    otherBranch.traverseAndDeformGeometries(false);

    outputTree = outputTree + outBranch;
    outputTree = outputTree + otherBranch;
    outputTree.resetInletOutlets();
}

double
SegmentationsMerger::
rotateBifurcation(Vector3D initialCenter,
                  Vector3D initialAxis,
                  Vector3D initialTransverse,
                  double initialRadius,
                  SegmentationParserPtr segmentationFather,
                  SegmentationParserPtr segmentationChild,
                  std::shared_ptr<BifurcationSymmetric> bifurcation)
{
    Vector3D refAxis = bifurcation->getInletNormal();
    Vector3D axis;
    double alpha;
    BuildingBlock::computeRotationAxisAndAngle(refAxis, initialAxis,
                                               axis, alpha);

    Matrix3D R = BuildingBlock::computeRotationMatrix(axis, alpha);
    Vector3D tr = initialCenter - R * bifurcation->getCenter();
    double scale = initialRadius / bifurcation->getInletRadius();

    Vector3D axisTransv = bifurcation->getTransverse();

    Vector3D axis_;
    double alphaTransv;
    BuildingBlock::computeRotationAxisAndAngle(R*axisTransv, initialTransverse,
                                               axis_, alphaTransv);

    int status;
    status = status || bifurcation->setParameterValue("bx", tr[0]);
    status = status || bifurcation->setParameterValue("by", tr[1]);
    status = status || bifurcation->setParameterValue("bz", tr[2]);
    status = status || bifurcation->setParameterValue("rotation_axis_x", axis[3]);
    status = status || bifurcation->setParameterValue("rotation_axis_y", axis[4]);
    status = status || bifurcation->setParameterValue("rotation_axis_z", axis[5]);
    status = status || bifurcation->setParameterValue("alpha", alpha);
    status = status || bifurcation->setParameterValue("scale", scale);

    std::vector<double> params(7, 0.0);
    params[0] = alphaTransv;

    return optimizeLoss(bifurcation, segmentationFather, segmentationChild,
                        params, 1e-3, 1e-7, 40000);
}

double
SegmentationsMerger::
findClosestPoint(const Vector3D& targetPoint,
                 SegmentationParserPtr toSearchIn,
                 unsigned int startIndexSearch,
                 unsigned int& indexOfClosestPoint)
{
    std::vector<Contour> contours = toSearchIn->getContours();
    unsigned int endIndex = toSearchIn->getIndexEnd();

    startIndexSearch = startIndexSearch < toSearchIn->getIndexBegin() ?
                       toSearchIn->getIndexBegin() : startIndexSearch;

    double minDist = 1e12;

    for (unsigned int i = startIndexSearch; i < endIndex; i++)
    {
        // here we rely on the fact that the distance should be monotone until
        // the minimum is reached
        bool found = false;
        double curDist = (contours[i].M_center-targetPoint).norm();
        if (curDist < minDist)
        {
            minDist = curDist;
            indexOfClosestPoint = i;
            found = true;
        }
        if (!found)
            break;
    }

    return minDist;
}

double
SegmentationsMerger::
placeBifurcation(Vector3D initialCenter,
                 Vector3D initialAxis,
                 Vector3D initialTransverse,
                 double initialRadius,
                 SegmentationParserPtr segmentationFather,
                 SegmentationParserPtr segmentationChild,
                 std::shared_ptr<BifurcationSymmetric> bifurcation)
{
    Vector3D refAxis = bifurcation->getInletNormal();
    Vector3D axis;
    double alpha;
    BuildingBlock::computeRotationAxisAndAngle(refAxis, initialAxis,
                                               axis, alpha);

    Matrix3D R = BuildingBlock::computeRotationMatrix(axis, alpha);
    Vector3D tr = initialCenter - R * bifurcation->getCenter();
    double scale = initialRadius / bifurcation->getInletRadius();

    Vector3D axisTransv = bifurcation->getTransverse();

    Vector3D axis_;
    double alphaTransv;
    BuildingBlock::computeRotationAxisAndAngle(R*axisTransv, initialTransverse,
                                               axis_, alphaTransv);

    std::vector<double> params(15, 0.0);
    params[0] = tr[0];
    params[1] = tr[1];
    params[2] = tr[2];
    params[3] = axis[0];
    params[4] = axis[1];
    params[5] = axis[2];
    params[6] = alpha;
    params[7] = scale;
    params[8] = -alphaTransv;

    return optimizeLoss(bifurcation, segmentationFather, segmentationChild,
                        params, 1e-5, 1e-7, 40000);
}

// params: bx,by,bz,rotation_axis_x,rotation_axis_y,rotation_axis_z,
// alpha, scale, alpha_axis, out1_alphax, out1_alphay, out1_alphaz,
// out2_alphax, out2_alphay, out2_alphaz
int
SegmentationsMerger::
deformBifurcation(std::shared_ptr<BifurcationSymmetric> bifurcation,
                  std::vector<double> params)
{
    int status = 0;

    auto map = bifurcation->getParametersMap();

    status = status || bifurcation->setParameterValue("bx", params[0]);
    status = status || bifurcation->setParameterValue("by", params[1]);
    status = status || bifurcation->setParameterValue("bz", params[2]);
    status = status || bifurcation->setParameterValue("rotation_axis_x", params[3]);
    status = status || bifurcation->setParameterValue("rotation_axis_y", params[4]);
    status = status || bifurcation->setParameterValue("rotation_axis_z", params[5]);
    status = status || bifurcation->setParameterValue("alpha", params[6]);
    status = status || bifurcation->setParameterValue("scale", params[7]);
    status = status || bifurcation->setParameterValue("alpha_axis", params[8]);
    status = status || bifurcation->setParameterValue("out1_alphax", params[9]);
    status = status || bifurcation->setParameterValue("out1_alphay", params[10]);
    status = status || bifurcation->setParameterValue("out1_alphaz", params[11]);
    status = status || bifurcation->setParameterValue("out2_alphax", params[12]);
    status = status || bifurcation->setParameterValue("out2_alphay", params[13]);
    status = status || bifurcation->setParameterValue("out2_alphaz", params[14]);

    bifurcation->applyGlobalTransformation(false);

    return status;
}

// params: alpha_axis, out1_alphax, out1_alphay, out1_alphaz,
// out2_alphax, out2_alphay, out2_alphaz
int
SegmentationsMerger::
deformPlacedBifurcation(std::shared_ptr<BifurcationSymmetric> bifurcation,
                        std::vector<double> params)
{
    int status = 0;

    auto map = bifurcation->getParametersMap();

    status = status || bifurcation->setParameterValue("alpha_axis", params[0]);
    status = status || bifurcation->setParameterValue("out1_alphax", params[1]);
    status = status || bifurcation->setParameterValue("out1_alphay", params[2]);
    status = status || bifurcation->setParameterValue("out1_alphaz", params[3]);
    status = status || bifurcation->setParameterValue("out2_alphax", params[4]);
    status = status || bifurcation->setParameterValue("out2_alphay", params[5]);
    status = status || bifurcation->setParameterValue("out2_alphaz", params[6]);

    bifurcation->applyGlobalTransformation(false);

    return status;
}

double
SegmentationsMerger::
optimizeLoss(std::shared_ptr<BifurcationSymmetric> bifurcation,
             SegmentationParserPtr segmentationFather,
             SegmentationParserPtr segmentationChild,
             std::vector<double>& params,
             double lambda,
             const double tol, const unsigned int nMax)
{
    printlog(MAGENTA, "[SegmentationsMerger] optimizing bifurcation position ...\n", M_verbose);

    const double eps = 1e-8;

    unsigned int nparams = params.size();
    std::vector<double> paramincr(nparams, 0.0);

    double incr = tol+1;
    unsigned int k = 0;

    double loss;
    int status;
    if (nparams == 7)
    {
        status = deformPlacedBifurcation(bifurcation, params);
        loss = computeLoss(bifurcation, segmentationFather, segmentationChild,0,1);
    }
    else
    {
        status = deformBifurcation(bifurcation, params);
        loss = computeLoss(bifurcation, segmentationFather, segmentationChild,1,1);
    }
    bifurcation->resetInletOutlets();
    double oldloss = loss;
    while (incr > tol && k <= nMax)
    {
        for (unsigned i = 0; i < nparams; i++)
        {
            params[i] = params[i] + eps;
            if (nparams == 7)
            {
                status = deformPlacedBifurcation(bifurcation, params);
                paramincr[i] = (computeLoss(bifurcation, segmentationFather,
                                           segmentationChild, 0, 1) - loss)/eps;
            }
            else
            {
                status = deformBifurcation(bifurcation, params);
                paramincr[i] = (computeLoss(bifurcation, segmentationFather,
                                           segmentationChild, 1, 1) - loss)/eps;
            }

            bifurcation->resetInletOutlets();

            params[i] = params[i] - eps;
        }

        for (unsigned int i = 0; i < nparams; i++)
        {
            params[i] -= lambda * paramincr[i];
        }

        if (nparams == 7)
        {
            status = deformPlacedBifurcation(bifurcation, params);
            loss = computeLoss(bifurcation, segmentationFather, segmentationChild,0,1);
        }
        else
        {
            status = deformBifurcation(bifurcation, params);
            loss = computeLoss(bifurcation, segmentationFather, segmentationChild,1,1);
        }
        bifurcation->resetInletOutlets();
        incr = std::abs((loss - oldloss)/loss);

        std::string msg;
        if (k % 1000 == 0)
        {
            msg = "it = " + std::to_string(k);
            std::ostringstream streamOb2;
            streamOb2 << loss;
            msg += " lf = " + streamOb2.str();
            streamOb2.str("");
            streamOb2.clear();
            streamOb2 << incr;
            msg += " increment = " + streamOb2.str();
            printlog(YELLOW, msg + "\n", M_verbose);
        }

        oldloss = loss;
        k++;
    }
    if (nparams == 7)
    {
        status = deformPlacedBifurcation(bifurcation, params);
        loss = computeLoss(bifurcation, segmentationFather, segmentationChild,0,1);
    }
    else
    {
        status = deformBifurcation(bifurcation, params);
        loss = computeLoss(bifurcation, segmentationFather, segmentationChild,1,1);
    }
    return loss;
}

double
SegmentationsMerger::
computeLoss(std::shared_ptr<BifurcationSymmetric> bifurcation,
            SegmentationParserPtr segmentationFather,
            SegmentationParserPtr segmentationChild,
            const double inletConst, const double outletConst)
{
    double loss = 0.0;

    double const1 = 1.0;
    double const2 = 1.0;
    double const3 = 2.0;

    Contour inlet = bifurcation->getInlet();
    Contour outlet1 = bifurcation->getOutlet(0);
    Contour outlet2 = bifurcation->getOutlet(1);

    unsigned int index;
    // inlet
    findClosestPoint(inlet.M_center, segmentationFather, 0, index);
    loss += const1 * (inlet.M_center - segmentationFather->getContour(index).M_center).norm();
    loss += const2 * std::abs(inlet.M_normal.dot(segmentationFather->getContour(index).M_normal)-1.0);
    loss += const3 * std::abs(inlet.M_radius - segmentationFather->getContour(index).M_radius);

    findClosestPoint(inlet.M_center, segmentationChild, 0, index);
    loss += const1 * (inlet.M_center - segmentationChild->getContour(index).M_center).norm();
    loss += const2 * std::abs(inlet.M_normal.dot(segmentationChild->getContour(index).M_normal)-1.0);
    loss += const3 * std::abs(inlet.M_radius - segmentationChild->getContour(index).M_radius);

    loss *= inletConst;

    unsigned int index1, index2;
    // outlet1
    double dist1 = findClosestPoint(outlet1.M_center, segmentationFather, 0, index1);
    double dist2 = findClosestPoint(outlet1.M_center, segmentationChild, 0, index2);

    const1 *= outletConst;
    const2 *= outletConst;

    if (dist1 < dist2)
    {
        loss += const1 * (outlet1.M_center - segmentationFather->getContour(index1).M_center).norm();
        loss += const2 * std::abs(outlet1.M_normal.dot(segmentationFather->getContour(index1).M_normal)-1.0);

        findClosestPoint(outlet2.M_center, segmentationChild, 0, index2);

        loss += const1 * (outlet2.M_center - segmentationChild->getContour(index2).M_center).norm();
        loss += const2 * std::abs(outlet2.M_normal.dot(segmentationChild->getContour(index2).M_normal)-1.0);
    }
    else
    {
        loss += const1 * (outlet1.M_center - segmentationChild->getContour(index2).M_center).norm();
        loss += const2 * std::abs(outlet1.M_normal.dot(segmentationChild->getContour(index2).M_normal)-1.0);

        findClosestPoint(outlet2.M_center, segmentationFather, 0, index1);

        loss += const1 * (outlet2.M_center - segmentationFather->getContour(index1).M_center).norm();
        loss += const2 * std::abs(outlet2.M_normal.dot(segmentationFather->getContour(index1).M_normal)-1.0);
    }

    return loss;
}

}  // namespace RedMA
