#include <ChebyshevBasisFunction.hpp>

namespace RedMA
{

ChebyshevBasisFunction::
ChebyshevBasisFunction(const GeometricFace& face,
                       unsigned int nMax) :
  BasisFunctionFunctor(face),
  M_sqrtPIm1(1.0/std::sqrt(M_PI))
{
    M_nMax = nMax;
    M_R = face.M_radius;

    for (int n = 0; n <= nMax; n++)
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
    unsigned int ind = (k * M_PI) / (n + 1);

    returnVal = M_sqrtPIm1 * chebyshevU(x/M_R*std::cos(ind) +
                                        y/M_R*std::sin(ind), n)/M_R;

    return returnVal;
}

}  // namespace RedMA
