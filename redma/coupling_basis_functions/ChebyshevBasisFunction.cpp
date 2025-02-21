#include "ChebyshevBasisFunction.hpp"

namespace RedMA
{

ChebyshevBasisFunction::
ChebyshevBasisFunction(const GeometricFace& face,
                       unsigned int nMax) :
  BasisFunctionFunctor(face),
  M_sqrtPIm1(1.0/std::sqrt(M_PI))
{
    M_nMax = nMax;

    for (int n = 0; n <= M_nMax; n++)
    {
        for (int k = 0; k <= n; k++)
        {
            M_ns.push_back(n);
            M_ks.push_back(k);
        }
    }

    M_nBasisFunctions = M_ns.size();
    M_type = "chebyshev";
}

double
ChebyshevBasisFunction::
chebyshevU(const double& x, const unsigned int& n)
{
    if (n == 0)
        return 1.0;
    else if (n == 1)
        return 2*x;
    return 2*x*chebyshevU(x,n-1) - chebyshevU(x,n-2);
}


ChebyshevBasisFunction::return_Type
ChebyshevBasisFunction::
operator()(const Vector3D& pos)
{
    double returnVal;

    double x;
    double y;
    getLocalXAndY(pos, x, y);

    unsigned int k = M_ks[M_index];
    unsigned int n = M_ns[M_index];
    double ind = static_cast<float>((k * M_PI)) / (n + 1);

    double R1 = M_face.M_radius1;
    double R2 = M_face.M_radius2;
    double R = std::sqrt(R1*R2);

    returnVal = M_sqrtPIm1 * chebyshevU(x/R1*std::cos(ind) +
                                        y/R2*std::sin(ind), n) / R;

    return returnVal;
}

}  // namespace RedMA
