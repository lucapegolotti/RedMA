#include "Tube.hpp"

namespace RedMA
{

Tube::
Tube(commPtr_Type comm, std::string refinement, bool verbose,
     int diameter, int length) :
  BuildingBlock(comm, refinement, verbose),
  M_length(length),
  M_diameter(diameter)
{
    M_name = "Tube";
    std::string ratio = std::to_string(diameter);
    ratio += "x";
    ratio += std::to_string(length);
    M_datafileName = "data_mesh";
    M_meshName = "tube_" + ratio + "/tube_" + ratio + "_" +
                 getStringMesh(refinement) + ".mesh";

    // center of inlet (reference configuration)
    M_inletCenterRef[0] = 0.0;
    M_inletCenterRef[1] = 0.0;
    M_inletCenterRef[2] = 0.0;

    // center of outlet (reference configuration)
    M_outletCenterRef[0] = 0.0;
    M_outletCenterRef[1] = 0.0;
    M_outletCenterRef[2] = length;

    // normal of inlet (reference configuration)
    M_inletNormalRef[0] = 0.0;
    M_inletNormalRef[1] = 0.0;
    M_inletNormalRef[2] = -1.0;

    // outlet of outlet (reference configuration)
    M_outletNormalRef[0] = 0.0;
    M_outletNormalRef[1] = 0.0;
    M_outletNormalRef[2] = 1.0;

    M_inletRadiusRef = (diameter * 1.0) / 2.0;
    M_outletRadiusRef = (diameter * 1.0) / 2.0;

    resetInletOutlets();

    // it is important to fill parametersHandler right at this level because then
    // the keys will be used in the parser to check the values in the XML file
    // center of inlet

    const bool randomizible = true;

    M_parametersHandler.registerParameter("bend", 0.0, 0, M_PI/2, randomizible, false);
    M_parametersHandler.registerParameter("L_ratio", 1.0, 0.7, 1.3, randomizible);
    M_parametersHandler.registerParameter("Rout_ratio", 1.0, 0.6, 1.4, randomizible);
    M_parametersHandler.registerParameter("use_linear_elasticity", 0.0, 0.0, 1.0);
}

void
Tube::
resetInletOutlets()
{
    GeometricFace inlet(M_inletCenterRef, M_inletNormalRef, M_inletRadiusRef, 1, 30);
    GeometricFace outlet(M_outletCenterRef, M_outletNormalRef, M_outletRadiusRef, 2, 31);

    M_inlet = inlet;
    M_outlets.clear();
    M_outlets.push_back(outlet);
}


std::string
Tube::
getOptionalParameter(unsigned int index)
{
    int retValue;
    // default diameter
    if (index == 0)
        retValue = static_cast<int>(M_inletRadiusRef * 2);
    // default length
    else if (index == 1)
        retValue = static_cast<int>(M_outletCenterRef[2]);
    else
        throw new Exception("Optional parameter does not exist!");
    return std::to_string(retValue);
}

std::string
Tube::
getStringMesh(std::string refinement)
{
    if (refinement[0] == 'h')
        return refinement;
    if (!std::strcmp(refinement.c_str(),"coarse"))
        return "h0.35";
    if (!std::strcmp(refinement.c_str(),"normal"))
        return "h0.25";
    if (!std::strcmp(refinement.c_str(),"fine"))
        return "h0.10";
    return "h0.25";
}


void
Tube::
applyNonAffineTransformation(bool transformMesh)
{
    printlog(MAGENTA, "[" + M_name +
                    " BuildingBlock] applying non affine transformation ...\n",
                    M_verbose);
    std::shared_ptr<LifeV::MeshUtility::MeshTransformer<mesh_Type> > transformer;

    if (transformMesh)
        transformer.reset(new LifeV::MeshUtility::MeshTransformer<mesh_Type>(*M_mesh));

    nonAffineScaling(M_parametersHandler["L_ratio"],
                     M_parametersHandler["Rout_ratio"],
                     transformer);

    bend(M_parametersHandler["bend"], transformer, transformMesh);
    printlog(MAGENTA, "done\n", M_verbose);
}

void
Tube::
nonAffineScaling(const double& lengthRatio, const double& radiusRatio,
                 std::shared_ptr<Transformer> transformer)
{
    std::string msg = std::string("[") + M_name + " BuildingBlock]";
    msg = msg + " scaling with lenghtRatio = " + std::to_string(lengthRatio) +
          " and radiusRatio = " + std::to_string(radiusRatio) + "\n";
    printlog(GREEN, msg, M_verbose);

    auto foo = std::bind(scalingFunction, std::placeholders::_1,
                         std::placeholders::_2, std::placeholders::_3,
                         lengthRatio, radiusRatio, M_outletCenterRef[2]);
    if (transformer)
        transformer->transformMesh(foo);

    M_outlets[0].M_radius = M_outlets[0].M_radius * radiusRatio;
    M_outlets[0].M_center[2] = M_outlets[0].M_center[2] * lengthRatio;
}

void
Tube::
scalingFunction(double& x, double& y, double& z,
                const double& lenghtRatio, const double& outRadiusRatio,
                const double& L)
{
    // 15 is the length of the tube
    double curRatio = 1. - (1. - outRadiusRatio) * z / L;
    z = z * lenghtRatio;
    x = x * curRatio;
    y = y * curRatio;
}

void
Tube::
bend(const double& bendAngle, std::shared_ptr<Transformer> transformer, bool transformMesh)
{
    std::string msg = std::string("[") + M_name + " BuildingBlock]";
    msg = msg + " bending with bendAngle = " + std::to_string(bendAngle) + "\n";
    printlog(GREEN, msg, M_verbose);

    if (bendAngle > 1e-5)
    {
        if (M_parametersHandler["use_linear_elasticity"])
        {
            using namespace std::placeholders;

            Vector3D rotationCenter = M_inlet.M_center;
            Matrix3D rotationMatrix = computeRotationMatrix(1, bendAngle);

            auto foo = std::bind(bendFunction,
                                 std::placeholders::_1,
                                 std::placeholders::_2,
                                 std::placeholders::_3,
                                 std::placeholders::_4,
                                 std::placeholders::_5, rotationCenter,
                                 rotationMatrix);

            if (transformMesh)
            {
                NonAffineDeformer nAffineDeformer(M_mesh, M_comm, M_verbose);

                LifeV::BCFunctionBase zeroFunction(BuildingBlock::fZero);
            	LifeV::BCFunctionBase outletFunction(foo);

                std::shared_ptr<LifeV::BCHandler> bcs(new LifeV::BCHandler);
                bcs->addBC("Inflow", 1, LifeV::Essential, LifeV::Full,
                           zeroFunction, 3);
                bcs->addBC("Outflow", 2, LifeV::Essential, LifeV::Full,
                           outletFunction, 3);

                CoutRedirecter ct;
                ct.redirect();
                nAffineDeformer.applyBCs(bcs);
                std::string xmlFilename = M_datafile("geometric_structure/xml_file",
                                                     "data/SolverParamList.xml");
                nAffineDeformer.setXMLsolver(xmlFilename);
                nAffineDeformer.deformMesh(*transformer);
                printlog(CYAN, ct.restore(), M_verbose);
            }

            // modify outlet
            Vector3D& center = M_outlets[0].M_center;
            Vector3D& normal = M_outlets[0].M_normal;

            center = rotationMatrix*(center - rotationCenter) + rotationCenter;
            normal = normal * rotationMatrix;
        }
        else
        {
            // this option gives smoother results
            // tube is laid in z direction (-> length == z component)
            auto foo = std::bind(bendFunctionAnalytic, std::placeholders::_1,
                                 std::placeholders::_2, std::placeholders::_3,
                                 bendAngle, M_outlets[0].M_center[2]);
            if (transformMesh)
                transformer->transformMesh(foo);
            // modify outlet
            Vector3D& center = M_outlets[0].M_center;
            Vector3D& normal = M_outlets[0].M_normal;

            bendFunctionAnalytic(center[0], center[1], center[2], bendAngle,
                                 M_outlets[0].M_center[2]);
            normal[0] = -std::sin(bendAngle);
            normal[1] = 0;
            normal[2] = std::cos(bendAngle);
        }
    }
}

double
Tube::
bendFunction(const double& t, const double& x, const double& y,
             const double& z, const LifeV::ID& i, const Vector3D& rotationCenter,
             const Matrix3D& rotationMatrix)
{
    Vector3D initial(x, y, z);
    Vector3D point(x, y, z);

    // bring to origin
    point = point - rotationCenter;

    // rotate point
    point = rotationMatrix * point;

    // bring back
    point = point + rotationCenter;

    return point[i] - initial[i];
}

void
Tube::
bendFunctionAnalytic(double &x, double &y, double &z,
                     const double& bendAngle, const double& L)
{
    // we know: lenght of the tube at initial configuration L and angle we want
    // alpha.
    // We imagine that the midline of the tube coincides with the arc of a circle.
    // Therefore, we have that
    // L = alpha * r => r = L / alpha, where R is the radius of the circle

    double r = L / bendAngle;

    // then, we determine the position of the current point in the arclenght
    double ratio = z / L;
    double actualAngle = bendAngle * ratio;

    // these refer to the point in the central line. Note that we assume that
    // the initial x of the central line is 0 because it shouldn't have been
    // moved up to this point
    double newX = -(r - r*std::cos(actualAngle));
    double newZ = r * std::sin(actualAngle);

    Vector3D xVersor;
    xVersor[0] = std::cos(actualAngle);
    xVersor[1] = 0;
    xVersor[2] = std::sin(actualAngle);

    newX = newX + x * xVersor[0];
    newZ = newZ + x * xVersor[2];
    x = newX;
    z = newZ;
}

}