#include "BifurcationSymmetric.hpp"

namespace RedMA
{

BifurcationSymmetric::
BifurcationSymmetric(EPETRACOMM comm, std::string refinement,
                     bool verbose, int angle, bool randomizable) :
  BuildingBlock(comm, refinement, verbose),
  M_angle(angle)
{
    M_name = "BifurcationSymmetric";

    M_datafileName = "bifurcation_symmetric_" + refinement + "_data";

    std::string alpha = std::to_string(angle);
    M_datafileName = "data_mesh";
    M_meshName = "bifurcation_symmetric/bif_sym_alpha" + alpha + "_" +
                 getStringMesh(refinement) + ".mesh";

    double angleRadiants = 2 * static_cast<double>(angle) / 360. * M_PI / 2;

    // values are referred to the reference configuration

    M_inletCenterRef[0] = 0;
    M_inletCenterRef[1] = 0;
    M_inletCenterRef[2] = 0;

    M_inletNormalRef[0] = 0;
    M_inletNormalRef[1] = 0;
    M_inletNormalRef[2] = -1;

    M_outlet1CenterRef[0] = 0;
    M_outlet1CenterRef[1] = 0.6;
    M_outlet1CenterRef[2] = 1.3;

    M_outlet1NormalRef[0] = 0.0;
    M_outlet1NormalRef[1] = std::sin(angleRadiants);
    M_outlet1NormalRef[2] = std::cos(angleRadiants);

    M_outlet2CenterRef[0] = 0;
    M_outlet2CenterRef[1] = -0.6;
    M_outlet2CenterRef[2] = 1.3;

    M_outlet2NormalRef[0] = 0;
    M_outlet2NormalRef[1] = -std::sin(angleRadiants);
    M_outlet2NormalRef[2] = std::cos(angleRadiants);

    M_inletRadiusRef = 0.5;
    M_outlet1RadiusRef = 0.5;
    M_outlet2RadiusRef = 0.5;

    M_transverse[0] = 1.0;
    M_transverse[1] = 0.0;
    M_transverse[2] = 0.0;

    resetInletOutlets();

    const double maxAngle = 0.4;
    M_parametersHandler.registerParameter("out1_alphax", 0.0, -maxAngle,
                                          maxAngle, randomizable);
    M_parametersHandler.registerParameter("out1_alphay", 0.0, -maxAngle,
                                          maxAngle, randomizable);
    M_parametersHandler.registerParameter("out1_alphaz", 0.0, -maxAngle,
                                          maxAngle, randomizable);
    M_parametersHandler.registerParameter("out2_alphax", 0.0, -maxAngle,
                                          maxAngle, randomizable);
    M_parametersHandler.registerParameter("out2_alphay", 0.0, -maxAngle,
                                          maxAngle, randomizable);
    M_parametersHandler.registerParameter("out2_alphaz", 0.0, -maxAngle,
                                          maxAngle, randomizable);
    computeCenter();

    M_identity3D(0,0) = 1;
    M_identity3D(0,1) = 0;
    M_identity3D(0,2) = 0;
    M_identity3D(1,0) = 0;
    M_identity3D(1,1) = 1;
    M_identity3D(1,2) = 0;
    M_identity3D(2,0) = 0;
    M_identity3D(2,1) = 0;
    M_identity3D(2,2) = 1;
}

std::string
BifurcationSymmetric::
getOptionalParameter(unsigned int index)
{
    int retValue;

    if (index == 0)
        retValue = static_cast<int>(M_angle);

    return std::to_string(retValue);
}


void
BifurcationSymmetric::
resetInletOutlets()
{
    GeometricFace inlet(M_inletCenterRef, M_inletNormalRef, M_inletRadiusRef, 1);
    GeometricFace outlet1(M_outlet1CenterRef, M_outlet1NormalRef, M_outlet1RadiusRef, 2);
    GeometricFace outlet2(M_outlet2CenterRef, M_outlet2NormalRef, M_outlet2RadiusRef, 3);

    M_inlet = inlet;

    M_outlets.clear();
    M_outlets.push_back(outlet1);
    M_outlets.push_back(outlet2);
}

void
BifurcationSymmetric::
computeCenter()
{
    double angleOutlet1 = std::atan(M_outlet1NormalRef[1] /
                                    M_outlet1NormalRef[2]);

    double zCenter = M_outlet1NormalRef[2] -
                     M_outlet1CenterRef[1] * std::tan(angleOutlet1);

    M_center[0] = 0;
    M_center[1] = 0;
    M_center[2] = zCenter;
}

std::string
BifurcationSymmetric::
getStringMesh(std::string refinement)
{
    return "h0.10";
}

void
BifurcationSymmetric::
bend(const double& out1_alphax, const double& out1_alphay, const double& out1_alphaz,
     const double& out2_alphax, const double& out2_alphay, const double& out2_alphaz,
     shp<Transformer> transformer, bool transformMesh)
{
    if ((std::abs(out1_alphax) > 0 || std::abs(out1_alphay) > 0 || std::abs(out1_alphaz) > 0) ||
        (std::abs(out2_alphax) > 0 || std::abs(out2_alphay) > 0 || std::abs(out2_alphaz) > 0))
    {
        std::string msg = std::string("[") + M_name + " BuildingBlock]";
        msg = msg + " bending with angles = (" + std::to_string(out1_alphax)
              + ", " + std::to_string(out1_alphay) +
                ", " + std::to_string(out1_alphaz) + ") at outlet1, and (" +
                       std::to_string(out2_alphax) + ", " +
                       std::to_string(out2_alphay) + ", " +
                       std::to_string(out2_alphaz) + ")"
              + " at outlet2" + "\n";
        printlog(GREEN, msg, M_verbose);

        using namespace std::placeholders;

        Vector3D rotationCenter = M_center;

        // handle rotation for outlet 1
        Matrix3D rotationMatrix =
                 computeRotationMatrix(0, out1_alphax) *
                 computeRotationMatrix(1, out1_alphay) *
                 computeRotationMatrix(2, out1_alphaz);

        Vector3D rotatedCenterOutlet1;
        Vector3D rotatedNormalOutlet1;
        rotateGeometricFace(M_outlets[0], rotatedCenterOutlet1, rotatedNormalOutlet1,
                            rotationMatrix, rotationCenter);

        auto fooOutlet1 = std::bind(outletMapFunction,
                                    std::placeholders::_1,
                                    std::placeholders::_2,
                                    std::placeholders::_3,
                                    std::placeholders::_4,
                                    std::placeholders::_5,
                                    M_outlets[0], rotatedCenterOutlet1,
                                    rotationMatrix);

        M_outlets[0].M_center = rotatedCenterOutlet1;
        M_outlets[0].M_normal = rotatedNormalOutlet1;

        // handle rotation for outlet 2
        rotationMatrix = computeRotationMatrix(0, out2_alphax) *
                         computeRotationMatrix(1, out2_alphay) *
                         computeRotationMatrix(2, out2_alphaz);

        Vector3D rotatedCenterOutlet2;
        Vector3D rotatedNormalOutlet2;
        rotateGeometricFace(M_outlets[1], rotatedCenterOutlet2, rotatedNormalOutlet2,
                            rotationMatrix, rotationCenter);

        auto fooOutlet2 = std::bind(outletMapFunction,
                                    std::placeholders::_1,
                                    std::placeholders::_2,
                                    std::placeholders::_3,
                                    std::placeholders::_4,
                                    std::placeholders::_5,
                                    M_outlets[1], rotatedCenterOutlet2,
                                    rotationMatrix);

        M_outlets[1].M_center = rotatedCenterOutlet2;
        M_outlets[1].M_normal = rotatedNormalOutlet2;

        if (transformMesh)
        {
            NonAffineDeformer nAffineDeformer(M_mesh, M_comm, M_verbose);

            LifeV::BCFunctionBase zeroFunction(BuildingBlock::fZero);
            LifeV::BCFunctionBase outletFunction1(fooOutlet1);
            LifeV::BCFunctionBase outletFunction2(fooOutlet2);

            shp<LifeV::BCHandler> bcs(new LifeV::BCHandler);
            bcs->addBC("Inflow", 1, LifeV::Essential, LifeV::Full,
                       zeroFunction, 3);
            bcs->addBC("Outflow", 2, LifeV::Essential, LifeV::Full,
                       outletFunction1, 3);
            bcs->addBC("Outflow", 3, LifeV::Essential, LifeV::Full,
                       outletFunction2, 3);

            CoutRedirecter ct;
            ct.redirect();
            nAffineDeformer.applyBCs(bcs);
            std::string xmlFilename = M_datafile("geometric_structure/xmldeformer",
                                                 "SolverParamList.xml");
            nAffineDeformer.setXMLsolver(xmlFilename);
            nAffineDeformer.deformMesh(*transformer);
            printlog(CYAN, ct.restore(), M_verbose);
        }
    }
}

double
BifurcationSymmetric::
outletMapFunction(const double& t, const double& x,
                  const double& y, const double& z,
                  const LifeV::ID& i, const GeometricFace& targetFace,
                  const Vector3D& desiredCenter, const Matrix3D& rotationMatrix)
{
    Vector3D cur(x, y, z);
    Vector3D diff = cur - targetFace.M_center;

    Vector3D modifiedDif = rotationMatrix * diff;
    Vector3D newPoint = desiredCenter + modifiedDif;

    return newPoint[i] - cur[i];
}

void
BifurcationSymmetric::
rotateGeometricFace(const GeometricFace& face, Vector3D& rotatedCenter,
                    Vector3D& rotatedNormal, const Matrix3D& rotationMatrix,
                    const Vector3D& rotationCenter)
{
    rotatedNormal = rotationMatrix * face.M_normal;
    rotatedCenter = rotationMatrix * (face.M_center - rotationCenter) +
                    rotationCenter;
}


void
BifurcationSymmetric::
applyNonAffineTransformation(bool transformMesh)
{
    printlog(MAGENTA, "[" + M_name +
                    " BuildingBlock] applying non affine transformation ...\n",
                    M_verbose);

    shp<LifeV::MeshUtility::MeshTransformer<MESH> > transformer;
    if (transformMesh)
        transformer.reset(new Transformer(*M_mesh));

    bend(M_parametersHandler["out1_alphax"],
         M_parametersHandler["out1_alphay"],
         M_parametersHandler["out1_alphaz"],
         M_parametersHandler["out2_alphax"],
         M_parametersHandler["out2_alphay"],
         M_parametersHandler["out2_alphaz"], transformer, transformMesh);

    printlog(MAGENTA, "done\n", M_verbose);
}

}
